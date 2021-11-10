from astropy.io import fits

class IO:

    def _to_fits(self, uniq, optional_kw_dict=None, pre_v2=False):
        """
        Serializes a MOC to the FITS format.

        Parameters
        ----------
        uniq : `numpy.ndarray`
            The array of HEALPix cells representing the MOC to serialize.
        optional_kw_dict : dict
            Optional keywords arguments added to the FITS header.
        pre_v2 : used only for ST-MOC FITS serialization (to ensure backward compatibility)

        Returns
        -------
        thdulist : `astropy.io.fits.HDUList`
            The list of HDU tables.
        """
        tbhdu = fits.BinTableHDU.from_columns(
            fits.ColDefs([
                fits.Column(name=self._fits_column_name, format=self._fits_format, array=uniq)
            ]))
        if pre_v2:
            tbhdu.header.update(self._fits_header_keywords_pre_v2)
        else:
            tbhdu.header.update(self._fits_header_keywords)

        if optional_kw_dict:
            for key in optional_kw_dict:
                tbhdu.header[key] = optional_kw_dict[key]

        thdulist = fits.HDUList([fits.PrimaryHDU(), tbhdu])
        return thdulist

    def serialize(self, format='fits', optional_kw_dict=None, pre_v2=False):
        """
        Serializes the MOC into a specific format.

        Possible formats are FITS, JSON and STRING

        Parameters
        ----------
        format : str
            'fits' by default. The other possible choice is 'json' or 'str'.
        optional_kw_dict : dict
            Optional keywords arguments added to the FITS header. Only used if ``format`` equals to 'fits'.

        Returns
        -------
        result : `astropy.io.fits.HDUList` or JSON dictionary
            The result of the serialization.
        """
        formats = ('fits', 'json', 'str')
        if format not in formats:
            raise ValueError('format should be one of %s' % (str(formats)))

        uniq = self._uniq_format()

        if format == 'fits':
            result = self._to_fits(uniq=uniq, optional_kw_dict=optional_kw_dict, pre_v2=pre_v2)
        elif format == 'str':
            result = self._to_str(uniq=uniq)
        else:
            # json format serialization
            result = self._to_json(uniq)
            # WARN: use the rust to_json
            # result = self._to_json(self._interval_set.nested)

        return result



    def write(self, path, format='fits', overwrite=False, optional_kw_dict=None, pre_v2=False):
        """
        Writes the MOC to a file.

        Format can be 'fits' or 'json', though only the fits format is officially supported by the IVOA.

        Parameters
        ----------
        path : str
            The path to the file to save the MOC in.
        format : str, optional
            The format in which the MOC will be serialized before being saved. Possible formats are "fits" or "json".
            By default, ``format`` is set to "fits".
        overwrite : bool, optional
            If the file already exists and you want to overwrite it, then set the  ``overwrite`` keyword. Default to False.
        optional_kw_dict : optional
            Optional keywords arguments added to the FITS header. Only used if ``format`` equals to 'fits'.
        """
        import warnings
        warnings.warn('This method is deprecated. Use MOC.save(path, "fits") instead!', DeprecationWarning)
        serialization = self.serialize(format=format, optional_kw_dict=optional_kw_dict, pre_v2=pre_v2)

        if format == 'fits':
            serialization.writeto(path, overwrite=overwrite)
        else:
            import os
            file_exists = os.path.isfile(path)

            if file_exists and not overwrite:
                raise OSError('File {} already exists! Set ``overwrite`` to '
                            'True if you want to replace it.'.format(path))

            if format == 'json':
                import json
                with open(path, 'w') as f_out:
                    f_out.write(
                        json.dumps(serialization, sort_keys=True, indent=2)
                    )
            elif format == 'str':
                with open(path, 'w') as f_out:
                    f_out.write(serialization)
