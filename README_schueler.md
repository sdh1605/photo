# photo.py fuer Schuelerinnen und Schueler

Dieses Skript hilft dir, aus zwei FITS-Bildern ein Sternfarben-Helligkeits-Diagramm (CMD) zu machen.

Kurz gesagt:
- Es sucht Sterne im Bild.
- Es misst ihre Helligkeit.
- Es vergleicht zwei Filter (meist B und V).
- Es zeichnet ein CMD.

## Was ist ein CMD?

Ein CMD zeigt:
- auf der x-Achse die Farbe: B-V
- auf der y-Achse die Helligkeit: V

Wichtig:
- kleine V-Werte bedeuten helle Sterne
- grosse V-Werte bedeuten schwache Sterne

## Schnellstart

Wenn du in GitHub Codespaces oder mit VS Code Dev Containers arbeitest, kannst du direkt die im Repository enthaltene Container-Umgebung nutzen.
Sie installiert automatisch die Python-Pakete aus `requirements.txt`.

```bash
python photo.py /pfad/zum/fits-ordner
```

## Haeufige Befehle

Nur bestimmte FITS-Dateien (0-basiert):

```bash
python photo.py /pfad/zum/fits-ordner --indices 0 1
```

Bilddarstellung als PNG speichern:

```bash
python photo.py /pfad/zum/fits-ordner --save frames.png
```

Nur Apertur:

```bash
python photo.py /pfad/zum/fits-ordner --photometry-method aperture
```

Apertur pro Stern anpassen:

```bash
python photo.py /pfad/zum/fits-ordner --photometry-method aperture --aperture-per-star
```

Nur PSF:

```bash
python photo.py /pfad/zum/fits-ordner --photometry-method psf
```

CMD speichern:

```bash
python photo.py /pfad/zum/fits-ordner --save-cm cmd.png
```

Apertur und PSF direkt vergleichen:

```bash
python photo.py /pfad/zum/fits-ordner --compare-cmd-methods --save-cmd-compare cmd_compare.png
```

Filter aus Kommandozeile erzwingen (Header ignorieren):

```bash
python photo.py /pfad/zum/fits-ordner --filters B V --force-filters
```

Sternprofile (hell/mittel/dunkel) plus PSF-Fit anzeigen:

```bash
python photo.py /pfad/zum/fits-ordner --plot-stellar-profiles --save-stellar-profiles profile.png
```

## Unterschied: Apertur vs. PSF (einfach erklaert)

## Apertur-Photometrie

Idee:
- Du legst einen Kreis um einen Stern.
- Alles Licht im Kreis wird addiert.
- Hintergrundlicht wird mit einem Ring ausserhalb abgeschaetzt und abgezogen.

Gut daran:
- einfach
- schnell
- funktioniert gut, wenn Sterne weit genug auseinander sind

Schwachstelle:
- wenn Sterne sehr nah beieinander liegen, kann Licht vermischt werden

## PSF-Photometrie

Idee:
- Statt nur einen Kreis zu summieren, passt das Programm ein Stern-Modell an.
- So kann es ueberlappende Sterne besser trennen.

Gut daran:
- besser in dichten Sternfeldern
- oft genauere Trennung naher Sterne

Schwachstelle:
- braucht mehr Rechenzeit
- haengt davon ab, wie gut das Modell passt

## Was zeigt der Vergleichsplot?

Wenn du --compare-cmd-methods nutzt:
- blaue Punkte: Apertur
- orange Punkte: PSF
- Linien: derselbe Stern in beiden Methoden

Zusaetzlich stehen im Plot/Terminal:
- RMS Delta(B-V): typische Farbabweichung
- RMS Delta(V): typische Helligkeitsabweichung
- Tendenz: ist PSF eher heller/dunkler und blauer/roeter?

## Kalibrierung (optional)

Mit --calibrate kannst du Referenzsterne auswaehlen.
Das Programm erlaubt jetzt auch mehrere Referenzsterne: Du klickst mehrere passende Sterne im Bild an und gibst fuer jeden Stern die B- und V-Helligkeit ein.
Danach werden die einzelnen Nullpunkte zusammengefasst und das CMD wird mit einer stabileren Kalibrierung berechnet.

Beispiel:

```bash
python photo.py /pfad/zum/fits-ordner --calibrate --calibrate-file calib.txt
```

## Dereddening (optional)

Mit E(B-V) kann die interstellare Roetung korrigiert werden.
Dann sind Farben und Helligkeiten physikalisch besser vergleichbar.

## Typische Probleme

photutils fehlt:

```bash
pip install photutils
```

Keine Matches gefunden:
- --maxsep etwas groesser machen
- --fwhm und --threshold anpassen
- Filterzuordnung pruefen (--filters oder --force-filters)

Zwei nahe helle Sterne werden als ein Stern erkannt:
- --detect-fwhm-scale kleiner setzen (z. B. 0.75)
- ggf. auch 0.65 testen
- wenn zu viele falsche Sterne erkannt werden: --threshold leicht erhoehen

Beispiel:

```bash
python photo.py /pfad/zum/fits-ordner --fwhm 3.0 --detect-fwhm-scale 0.75
```

Zu viele Farbausreisser im CMD:
- --outlier-sigma groesser setzen (z.B. 4.0 oder 5.0)
