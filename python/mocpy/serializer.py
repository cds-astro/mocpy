from astropy.io import fits
from . import mocpy


class IO:
    def serialize(self, format="fits", optional_kw_dict=None, pre_v2=False):
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
        formats = ("fits", "json", "str")
        if format not in formats:
            raise ValueError("format should be one of %s" % (str(formats)))

        if format == "fits":
            hdulist = fits.HDUList.fromstring(
                mocpy.to_fits_raw(self._store_index, pre_v2)
            )
            hdu = hdulist[1]
            if optional_kw_dict:
                for key in optional_kw_dict:
                    hdu.header[key] = optional_kw_dict[key]
            return hdulist

        elif format == "str":
            result = self.to_string(format="ascii", fold=0)
        else:
            import json

            json_str = self.to_string(format="json")
            result = json.loads(json_str)

        return result

    def write(
        self, path, format="fits", overwrite=False, optional_kw_dict=None, pre_v2=False
    ):
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

        warnings.warn(
            'This method is deprecated. Use MOC.save(path, "fits") instead!',
            DeprecationWarning,
        )
        serialization = self.serialize(
            format=format, optional_kw_dict=optional_kw_dict, pre_v2=pre_v2
        )

        if format == "fits":
            serialization.writeto(path, overwrite=overwrite)
        else:
            import os

            file_exists = os.path.isfile(path)

            if file_exists and not overwrite:
                raise OSError(
                    "File {} already exists! Set ``overwrite`` to "
                    "True if you want to replace it.".format(path)
                )

            if format == "json":
                import json

                with open(path, "w") as f_out:
                    f_out.write(json.dumps(serialization, sort_keys=True, indent=2))
            elif format == "str":
                with open(path, "w") as f_out:
                    f_out.write(serialization)
