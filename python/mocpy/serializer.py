import warnings

from astropy.io import fits

from . import mocpy


class IO:
    """Input and outputs for MOCs."""

    def serialize(
        self,
        format="fits",
        optional_kw_dict=None,
        v1_compatible=False,
        pre_v2=None,
    ):
        """
        Serialize the MOC into a specific format.

        Possible formats are FITS, JSON and STRING

        Parameters
        ----------
        format : str
            'fits' by default. The other possible choice is 'json' or 'str'.
        optional_kw_dict : dict
            Optional keywords arguments added to the FITS header. Only used if ``format`` equals to 'fits'.
        v1_compatible : bool, optional
            If True, write a v2 MOC that can be read by v1 parsers (using NUNIQ ordering for spatial MOCs).
            Default to False (MOC 2.0 standard range-encoded).
        pre_v2 : bool, optional
            Deprecated since v0.21.0. Use ``v1_compatible`` instead.

        Returns
        -------
        result : `astropy.io.fits.HDUList` or JSON dictionary
            The result of the serialization.
        """
        formats = ("fits", "json", "str")
        if format not in formats:
            raise ValueError(f"format should be one of {formats}")

        if pre_v2 is not None:
            warnings.warn(
                "The 'pre_v2' argument is deprecated. Use 'v1_compatible' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            v1_compatible = pre_v2

        if format != "fits" and v1_compatible:
            warnings.warn(
                "The 'v1_compatible' argument is only valid when 'fits' is the output format. It will be ignored.",
                UserWarning,
                stacklevel=2,
            )

        if format == "fits":
            hdulist = fits.HDUList.fromstring(
                mocpy.to_fits_raw(self.store_index, v1_compatible),
            )
            hdu = hdulist[1]
            if optional_kw_dict:
                for key in optional_kw_dict:
                    hdu.header[key] = optional_kw_dict[key]
            return hdulist

        if format == "str":
            return self.to_string(format="ascii", fold=0)

        import json

        json_str = self.to_string(format="json")
        return json.loads(json_str)

    def write(
        self,
        path,
        format="fits",
        overwrite=False,
        optional_kw_dict=None,
        v1_compatible=False,
        pre_v2=None,
    ):
        """
        Write the MOC to a file.

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
        v1_compatible : bool, optional
            If True, write a v2 MOC that can be read by v1 parsers (using NUNIQ ordering for spatial MOCs).
            Default to False (MOC 2.0 standard range-encoded).
        pre_v2 : bool, optional
            Deprecated since v0.21.0. Use ``v1_compatible`` instead.
        """
        warnings.warn(
            'This method is deprecated. Use MOC.save(path, "fits") instead!',
            DeprecationWarning,
            stacklevel=2,
        )
        self.save(
            path,
            format,
            overwrite,
            v1_compatible=v1_compatible,
            fits_keywords=optional_kw_dict,
            pre_v2=pre_v2,
        )
