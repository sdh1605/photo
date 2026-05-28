# photo.py Cheat Sheet (10 Zeilen)
1. Start: python photo.py /pfad/zum/ordner
2. Apertur: python photo.py /pfad --photometry-method aperture
3. PSF: python photo.py /pfad --photometry-method psf
4. CMD speichern: python photo.py /pfad --save-cm cmd.png
5. Methodenvergleich: python photo.py /pfad --compare-cmd-methods
6. Vergleich speichern: python photo.py /pfad --compare-cmd-methods --save-cmd-compare cmp.png
7. Kalibrieren: python photo.py /pfad --calibrate --calibrate-file calib.txt
8. Deredden/Isochrone: python photo.py /pfad --calibrate --fit-isochrone --ebv 0.05
9. Problem keine Matches: --maxsep erhöhen, --fwhm/--threshold anpassen
10. Merke: Apertur = Kreis-Summe, PSF = Modellfit (besser bei dichten Sternfeldern)
