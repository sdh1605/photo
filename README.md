# photo.py

CLI-Tool zur Photometrie von zwei FITS-Bildern mit anschliessendem Color-Magnitude-Diagramm (CMD).

Das Skript:
- findet zwei FITS-Dateien in einem Ordner,
- detektiert Sterne,
- fuehrt Photometrie durch (Apertur oder PSF),
- matched Sterne zwischen zwei Filtern (typisch B und V),
- erstellt CMD-Plots,
- kann kalibrieren, deredden und Isochronen anpassen,
- kann Apertur vs. PSF direkt vergleichen.

## Voraussetzungen

Benötigte Python-Pakete:
- numpy
- matplotlib
- astropy
- photutils

Beispiel-Installation:

```bash
pip install numpy matplotlib astropy photutils
```

## Schnellstart

```bash
python photo.py /pfad/zum/fits-ordner
```

Mit expliziter Methode:

```bash
python photo.py /pfad/zum/fits-ordner --photometry-method aperture
python photo.py /pfad/zum/fits-ordner --photometry-method psf
```

CMD speichern:

```bash
python photo.py /pfad/zum/fits-ordner --save-cm cmd.png
```

Vergleich Apertur vs. PSF als gemeinsamer CMD-Plot:

```bash
python photo.py /pfad/zum/fits-ordner --compare-cmd-methods --save-cmd-compare cmd_compare.png
```

## Wichtige Optionen

- `--photometry-method {aperture,psf}`
  - Waehlt die Photometrie-Methode.
- `--aperture FLOAT`
  - Aperturradius in Pixeln. Ohne Angabe wird er automatisch aus FWHM berechnet.
- `--fwhm FLOAT`
  - FWHM fuer Quellenfindung.
- `--threshold FLOAT`
  - Sigma-Schwelle fuer Quellenfindung.
- `--maxsep FLOAT`
  - Maximaler Matching-Abstand in Pixeln.
- `--filters B V`
  - Filterzuordnung, falls keine gute Header-Info vorliegt.
- `--force-filters`
  - Erzwingt die `--filters`-Zuordnung.
- `--cluster-radius`, `--cluster-center`, `--annulus-radii`
  - Raeumliche Selektion von Sternen.
- `--calibrate`, `--calibrate-file`
  - Referenzstern-Kalibrierung laden/speichern.
- `--fit-isochrone`, `--ebv`, `--age-range`, `--age`
  - Isochronen-Anpassung fuer Entfernung/Alter.
- `--compare-cmd-methods`, `--save-cmd-compare`
  - Direkter Methodenvergleich im CMD.

## Workflow (kurz)

1. Zwei FITS-Dateien laden.
2. Sterne detektieren.
3. Photometrie pro Bild berechnen.
4. Sterne zwischen beiden Bildern robust matchen.
5. B und V zuordnen.
6. CMD erstellen (instrumentell oder kalibriert).
7. Optional: Dereddening und Isochrone.
8. Optional: Methodenvergleich Apertur vs. PSF mit Linien und Abweichungs-Statistik.

## Unterschied: Apertur vs. PSF-Photometrie

### Apertur-Photometrie

Prinzip:
- Summiert den Sternfluss in einem Kreis um die Sternposition.
- Schaetzt den lokalen Hintergrund in einem Ring (Annulus).
- Subtrahiert Hintergrundanteil von der Apertur-Summe.

Vorteile:
- Einfach, schnell, robust.
- Gut fuer eher isolierte Sterne und moderate Felddichte.

Nachteile:
- Bei ueberlappenden Sternen wird Fluss leicht vermischt.
- Abhaengiger von Radius-Wahl und Crowding.

### PSF-Photometrie

Prinzip:
- Passt ein PSF-Modell (hier Gaussian-basiert) an Sternprofile an.
- Trennt ueberlappende Quellen besser durch Modellanpassung.
- Liefert modellbasierten Sternfluss (`flux_fit`).

Vorteile:
- Besser bei dichtem Sternfeld (Crowding).
- Bessere Trennung naher Sterne.

Nachteile:
- Rechenaufwendiger.
- Abhaengig von der Guete des PSF-Modells.
- Bei Problemen wird im Skript auf Apertur-Photometrie zurueckgefallen.

## Methodenvergleich im Skript

Mit `--compare-cmd-methods` wird ein gemeinsamer CMD-Plot erzeugt:
- Apertur- und PSF-Punkte in unterschiedlichen Farben.
- Verbindungssegmente zwischen identischen Sternen.
- Kennzahlen aus den Differenzen (PSF minus Apertur):
  - RMS Delta(B-V)
  - RMS Delta(V)
  - mittlere Tendenz in Farbe (blauer/roeter)
  - mittlere Tendenz in Helligkeit (heller/dunkler)

Falls Kalibrierung aktiv ist, wird der Vergleich ebenfalls kalibriert berechnet.
Falls `E(B-V)` verfuegbar ist, erfolgt zusaetzlich dereddening im Vergleichsplot.

## Ausgabedateien

Je nach Optionen werden u.a. erzeugt:
- Bildplot der Frames
- CMD-Plot
- CMD-Methodenvergleich
- `calib.txt` (Kalibrierungsdaten)
- `photometry.txt` (Sternliste mit V und B-V)
- `cluster_parameters.txt` (bei Isochronen-Fit)
- `isochrone_fit.png` (bei Isochronen-Fit)

Standardmaessig werden relative Ausgabepfade im FITS-Ordner gespeichert.

## Typische Probleme

- `photutils` nicht installiert:
  - Fehlermeldung bei Import oder beim Start der Photometrie.
  - Loesung: `pip install photutils`.
- Keine Matches gefunden:
  - `--maxsep` vergroessern,
  - `--fwhm` und `--threshold` an Daten anpassen,
  - Filterzuordnung pruefen (`--filters`, `--force-filters`).
- Unerwartete Farben im CMD:
  - Kalibrierung pruefen,
  - B/V-Zuordnung pruefen,
  - outlier-Schwelle mit `--outlier-sigma` anpassen.
