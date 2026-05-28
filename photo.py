#!/usr/bin/env python3
"""
photo.py

CLI-Tool (Deutsch): Suche in einem Verzeichnis nach zwei FITS-Dateien, lade sie und zeige sie an.

Verwendung:
  python photo.py /pfad/zum/ordner

Optionen:
  -i, --indices   Indices (0-basiert) der anzuzeigenden FITS-Dateien
  -s, --save      Optionaler Pfad zum Speichern der Abbildung (PNG)
  -v, --verbose   Mehr Ausgaben
  --photometry-method METHOD  Photometrie-Methode: 'aperture' oder 'psf' (Standard: aperture)
  --aperture FLOAT    Aperturradius in Pixel (wenn nicht angegeben: automatisch 3.0×FWHM, empfohlen!)
  --fwhm FLOAT       FWHM für Quellenfindung (Standard: 3.0)
  --threshold FLOAT  Schwelle in Sigma für Quellenfindung (Standard: 5.0)
  --maxsep FLOAT     Max. Abstand beim Matching in Pixel (Standard: 5.0)
  --filters B V      Labels der Filter (Standard: B V)
  --save-cm FILE     Speicherpfad für das Color-Magnitude-Diagramm (relativ zum FITS-Ordner oder absolut)
    --compare-cmd-methods  Vergleicht CMDs aus aperture- und psf-Auswertung in einem gemeinsamen Plot
    --save-cmd-compare FILE  Speicherpfad für den CMD-Methodenvergleich (relativ zum FITS-Ordner oder absolut)
  --cluster-radius FLOAT  Radius in Pixel um das Zentrum für Sternhaufen-Filterung (Standard: 40% der Bildhöhe, 0 = deaktiviert)
  --cluster-center X Y    Zentrum (x, y) in Pixel für Sternhaufen-Filterung (Standard: Bildmitte)
  --annulus-radii R_INNER R_OUTER  Nur Sterne in einem Ring (Annulus) zwischen innerem (R_INNER) und
                          äußerem (R_OUTER) Radius in Pixel um das Zentrum berücksichtigen.
                          Überschreibt --cluster-radius wenn angegeben.
  --calibrate        Aktiviere interaktive Kalibrierung mit Referenzstern (klicken + Magnitudeneingabe)
  --calibrate-file FILE   Pfad zur Kalibrierungsdatei (relativ zum FITS-Ordner oder absolut).
                          Mit --calibrate: Speichern. Ohne --calibrate: Laden. Standard: calib.txt im FITS-Ordner.
  --fit-isochrone    Bestimme Entfernung und Alter durch Isochronen-Anpassung (benötigt --calibrate)
  --ebv FLOAT        Rötung E(B-V) in mag für Isochronen-Anpassung (wird interaktiv abgefragt falls nicht angegeben)
  --age-range MIN MAX    Altersbereich als log(Alter/Jahre) (Standard: 6.5 8.5 = 3-300 Millionen Jahre)
  --age FLOAT        Fixiertes Alter in Millionen Jahren (z.B. 10.6). Wenn gesetzt, wird nur die Entfernung bestimmt.

Das Skript sucht nach Dateien mit Endungen .fits oder .fit und versucht, die erste 2D-Image-HDU
anzuzeigen. Falls ein HDU 3D-Daten enthält, wird die erste Ebene verwendet.

Alle Ausgabedateien (Plots, Kalibrierungsdaten) werden standardmäßig im selben Verzeichnis wie
die FITS-Dateien gespeichert, außer es wird ein absoluter Pfad angegeben.

Isochronen-Anpassung:
  Das Programm lädt Isochrone von CMD 3.7 (http://stev.oapd.inaf.it) und bestimmt durch Chi²-Minimierung
  das beste Alter und die Entfernung des Sternhaufens. Die Rötung E(B-V) muss bekannt sein.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
import urllib.request
import urllib.parse
from typing import List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import ImageNormalize, ZScaleInterval

# photutils wird für Photometrie verwendet (separates Paket, nicht Teil von astropy)
try:
    from photutils.detection import DAOStarFinder
    from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
    from astropy.stats import sigma_clipped_stats
    # PSF-Photometrie Imports (optional)
    try:
        from photutils.psf import IterativePSFPhotometry
        HAVE_PSF = True
    except ImportError:
        try:
            from photutils.psf import PSFPhotometry
            HAVE_PSF = True
        except ImportError:
            HAVE_PSF = False
    HAVE_PHOTUTILS = True
except ImportError as e:
    HAVE_PHOTUTILS = False
    HAVE_PSF = False
    PHOTUTILS_ERROR = str(e)



def find_fits_files(folder: str) -> List[str]:
	"""Finde alle .fits oder .fit Dateien im Verzeichnis und gebe sie als sortierte Liste zurück."""
	all_entries = os.listdir(folder)
	fits_files = []
	for entry in all_entries:
		if entry.lower().endswith('.fits') or entry.lower().endswith('.fit'):
			fits_files.append(os.path.join(folder, entry))
	fits_files.sort()
	return fits_files


def get_output_path(fits_folder: str, filename: str) -> str:
	"""Generiere Ausgabepfad im selben Verzeichnis wie die FITS-Dateien.
	
	Args:
		fits_folder: Verzeichnis mit den FITS-Dateien
		filename: Name der Ausgabedatei (z.B. 'calib.txt', 'cmd.png')
	
	Returns:
		Vollständiger Pfad zur Ausgabedatei
	"""
	return os.path.join(fits_folder, filename)


def _get_image_data_from_hdul(hdul) -> Optional[np.ndarray]:
	"""Finde die erste HDU mit 2D-Image-Daten. Bei 3D: erste Scheibe zurückgeben."""
	for h in hdul:
		data = h.data
		if data is None:
			continue
		arr = np.asarray(data)
		if arr.ndim == 2:
			return arr
		if arr.ndim == 3 and arr.shape[0] > 0:
			return arr[0]
	return None


def assign_filters_by_flux(tab1, tab2):
    """Fallback: Bestimme Filter-Zuordnung basierend auf mittlerem Flux."""
    mean_flux_1 = np.mean(tab1['aperture_sum_bkgsub']) if tab1 is not None else 0
    mean_flux_2 = np.mean(tab2['aperture_sum_bkgsub']) if tab2 is not None else 0
    
    print(f"\nFlux-Analyse:")
    print(f"  Bild 1: mittlerer Flux = {mean_flux_1:.1f}")
    print(f"  Bild 2: mittlerer Flux = {mean_flux_2:.1f}")
    
    # V-Filter sollte typischerweise höheren Flux zeigen
    if mean_flux_1 > mean_flux_2 * 1.1:
        print("  → Bild 1 ist wahrscheinlich V (heller)")
        print("  → Zuordnung: tabB = Bild 2, tabV = Bild 1")
        return tab2, tab1  # (tabB, tabV)
    elif mean_flux_2 > mean_flux_1 * 1.1:
        print("  → Bild 2 ist wahrscheinlich V (heller)")
        print("  → Zuordnung: tabB = Bild 1, tabV = Bild 2")
        return tab1, tab2  # (tabB, tabV)
    else:
        print("  → Flüsse ähnlich, verwende Bild 1 = B, Bild 2 = V")
        return tab1, tab2  # (tabB, tabV)


def get_filter_from_header(header):
    """Extrahiere Filter-Information aus FITS-Header.
    
    Sucht nach gängigen Header-Keywords für Filter.
    Returns: Filter-Name (z.B. 'B', 'V') oder None
    """
    if header is None:
        print("    DEBUG: Header ist None")
        return None
    
    # Gängige FITS-Keywords für Filter (in Prioritätsreihenfolge)
    filter_keywords = ['FILTER', 'FILTER1', 'FILTNAM', 'FILT', 'INSFILTE', 'HIERARCH ESO INS FILT1 NAME']
    
    print(f"    DEBUG: Durchsuche {len(header)} Header-Einträge nach Filter-Keywords...")
    for keyword in filter_keywords:
        if keyword in header:
            filter_value = str(header[keyword]).strip().upper()
            print(f"    DEBUG: Gefunden {keyword} = '{filter_value}'")
            # Extrahiere nur den Buchstaben (z.B. "V" aus "V-Band" oder "SDSS-V")
            if filter_value:
                # Suche nach B, V, R, I, U, etc.
                for char in filter_value:
                    if char in 'BVRIUJHK':
                        print(f"    DEBUG: Filter '{char}' extrahiert aus '{filter_value}'")
                        return char
                # Falls kein einzelner Buchstabe gefunden, gebe den ganzen String zurück
                print(f"    DEBUG: Kein Standard-Filter gefunden, gebe '{filter_value}' zurück")
                return filter_value
    
    print(f"    DEBUG: Kein Filter-Keyword gefunden. Zeige erste 20 Header-Einträge:")
    for i, key in enumerate(sorted(header.keys())):
        if i >= 20:
            break
        print(f"      {key} = {header[key]}")
    
    return None


def load_fits_image(path: str) -> Tuple[Optional[np.ndarray], Optional[dict]]:
	"""Öffne eine FITS-Datei und liefere (data, header) zurück. data kann None sein."""
	try:
		with fits.open(path, memmap=False) as hdul:
			data = _get_image_data_from_hdul(hdul)
			# Behalte den vollen Header als dict
			header = dict(hdul[0].header) if len(hdul) > 0 and hdul[0].header is not None else None
			return data, header
	except Exception as e:
		raise RuntimeError(f"Fehler beim Öffnen von {path}: {e}")


def plot_images(images: List[Tuple[Optional[np.ndarray], Optional[dict]]], paths: List[str], 
              phot_tables: Optional[List] = None, aperture_radius: Optional[float] = None,
              cluster_radius: Optional[float] = None, cluster_center: Optional[Tuple[float, float]] = None,
              image_offset: Optional[np.ndarray] = None,
              matched_indices: Optional[List[Tuple[int, int]]] = None,
              annulus_inner: Optional[float] = None,
              save: Optional[str] = None) -> None:
    """Zeige FITS-Bilder mit Aperturen und optionalem Cluster-Radius oder Annulus.
    
    Args:
        image_offset: Offset [dx, dy] zwischen erstem und zweitem Bild (für Cluster-Kreis Verschiebung)
        matched_indices: Liste von (index1, index2) Tupeln für gematchte Sterne
        annulus_inner: Innerer Annulus-Radius in Pixel (wenn gesetzt, wird cluster_radius als äußerer Radius verwendet)
    """
    n = len(images)
    if n == 0:
        print("Keine Bilder zum Anzeigen.")
        return
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 6))
    if n == 1:
        axes = [axes]

    for i, (ax, (data, header), path) in enumerate(zip(axes, images, paths)):
        ax.set_title(os.path.basename(path))
        if data is None:
            ax.text(0.5, 0.5, 'Keine 2D-Daten gefunden', ha='center', va='center')
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        try:
            norm = ImageNormalize(data, interval=ZScaleInterval())
        except Exception:
            norm = None
        im = ax.imshow(data, cmap='gray', origin='lower', norm=norm)
        ax.set_xlabel(f"Shape: {data.shape}")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        
        # Bestimme das Zentrum für dieses spezifische Bild
        # Erstes Bild: Original-Zentrum
        # Zweites Bild: Zentrum + Offset
        current_center = cluster_center
        if cluster_center is not None and i > 0 and image_offset is not None:
            # Verschiebe das Zentrum für das zweite Bild
            current_center = (cluster_center[0] + image_offset[0], 
                            cluster_center[1] + image_offset[1])
        
        # Zeige Aperturen falls verfügbar
        if phot_tables and i < len(phot_tables) and phot_tables[i] is not None and aperture_radius:
            tab = phot_tables[i]
            
            # Wenn matched_indices verfügbar ist, zeige NUR gematchte Sterne
            if matched_indices is not None and len(matched_indices) > 0:
                if i == 0:
                    # Erstes Bild: zeige Sterne mit Indizes aus matched_indices[:,0]
                    indices_to_show = [idx1 for idx1, idx2 in matched_indices]
                else:
                    # Zweites Bild: zeige Sterne mit Indizes aus matched_indices[:,1]
                    indices_to_show = [idx2 for idx1, idx2 in matched_indices]
                
                coords = np.vstack((np.array(tab['x']), np.array(tab['y']))).T
                coords_to_show = coords[indices_to_show]
                positions = list(zip(coords_to_show[:, 0], coords_to_show[:, 1]))
            else:
                # Fallback: zeige alle Sterne (wie vorher)
                coords = np.vstack((np.array(tab['x']), np.array(tab['y']))).T
                
                # Filtere nach Radius/Annulus falls angegeben (verwende aktuelles Zentrum)
                if cluster_radius is not None and current_center is not None:
                    dist = np.linalg.norm(coords - np.array(current_center), axis=1)
                    if annulus_inner is not None:
                        mask = (dist >= annulus_inner) & (dist <= cluster_radius)
                    else:
                        mask = dist <= cluster_radius
                    coords_filtered = coords[mask]
                    positions = list(zip(coords_filtered[:, 0], coords_filtered[:, 1]))
                else:
                    positions = list(zip(tab['x'], tab['y']))
            
            apertures = CircularAperture(positions, r=aperture_radius)
            apertures.plot(color='red', lw=1.0, alpha=0.7, ax=ax)
        
        # Zeige Cluster-Radius-Kreis oder Annulus falls definiert (verwende aktuelles Zentrum)
        if cluster_radius is not None and current_center is not None:
            from matplotlib.patches import Circle
            # Äußerer Kreis
            circle = Circle(current_center, cluster_radius, fill=False, 
                          edgecolor='cyan', linewidth=2, linestyle='--', alpha=0.8)
            ax.add_patch(circle)
            # Innerer Kreis (nur bei Annulus)
            if annulus_inner is not None and annulus_inner > 0:
                inner_circle = Circle(current_center, annulus_inner, fill=False,
                                      edgecolor='cyan', linewidth=2, linestyle=':', alpha=0.8)
                ax.add_patch(inner_circle)
            # Markiere das Zentrum
            ax.plot(current_center[0], current_center[1], 'c+', markersize=15, markeredgewidth=2)

    plt.tight_layout()
    if save:
        try:
            fig.savefig(save, dpi=150)
            print(f"Abbildung gespeichert: {save}")
        except Exception as e:
            print(f"Fehler beim Speichern: {e}", file=sys.stderr)
    plt.show()


def detect_and_photometry(data: np.ndarray, fwhm: float = 3.0, threshold_sigma: float = 5.0, aperture_radius: float = None):
    """Erkenne Quellen mit DAOStarFinder und führe Aperturphotometrie mit Hintergrundsubtraktion aus.

    Liefert (phot_table, aperture_radius):
    - phot_table: Tabelle mit x, y, aperture_sum (Roh-Flux), sky_median (Himmelshintergrund pro Pixel),
                  und aperture_sum_bkgsub (hintergrundkorrigierter Flux).
    - aperture_radius: Der tatsächlich verwendete Aperturradius (wichtig wenn automatisch berechnet)
    
    Die Hintergrundsubtraktion erfolgt mittels Sky-Annulus (Kreisring um jede Quelle).
    
    Parameters:
        aperture_radius: Wenn None, wird automatisch 3.0×FWHM verwendet (empfohlen!).
                        Manueller Wert nur wenn spezielle Anforderungen vorliegen.
    
    Returns:
        (phot_table, aperture_radius): Tuple aus Photometrie-Tabelle und verwendetem Aperturradius
    """
    if not HAVE_PHOTUTILS:
        msg = (
            "Das Paket 'photutils' ist nicht installiert. Bitte installieren Sie es mit:\n"
            "pip install photutils\n"
            f"Original Fehler: {PHOTUTILS_ERROR}"
        )
        raise RuntimeError(msg)
    
    # Automatische Apertur-Berechnung basierend auf FWHM
    # Verwende 3.0×FWHM um auch helle (gesättigte) Sterne vollständig zu erfassen
    # Helle Sterne haben oft größere PSF durch Sättigung, Blooming, etc.
    # Für einheitliche Photometrie: eine Apertur für alle Sterne (basierend auf hellsten Sternen)
    if aperture_radius is None:
        aperture_radius = 3.0 * fwhm  # Konservativ: erfasst ~95% des Lichts
        print(f"  Automatische Apertur: radius = 3.0 × FWHM = 3.0 × {fwhm:.1f} = {aperture_radius:.1f} Pixel")
        print(f"    (3.0×FWHM ist konservativ und erfasst auch helle/gesättigte Sterne vollständig)")
    else:
        print(f"  Manuelle Apertur: radius = {aperture_radius:.1f} Pixel")
    
    # Robust stats für Quellenfindung
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    threshold = median + threshold_sigma * std

    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)
    sources = daofind(data - median)
    if sources is None or len(sources) == 0:
        return None

    # Verwende Zentren der Quellen für Aperturphotometrie
    positions = list(zip(sources['xcentroid'], sources['ycentroid']))
    
    # Signal-Apertur
    apertures = CircularAperture(positions, r=aperture_radius)
    phot_table = aperture_photometry(data, apertures)
    
    # Hintergrund-Annulus (typischerweise zwischen r_in = 1.5*r und r_out = 2.5*r)
    r_in = aperture_radius * 1.5
    r_out = aperture_radius * 2.5
    annulus_apertures = CircularAnnulus(positions, r_in=r_in, r_out=r_out)
    annulus_masks = annulus_apertures.to_mask(method='center')
    
    # Berechne den Median des Himmels in jedem Annulus
    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(data)
        annulus_data_1d = annulus_data[mask.data > 0]
        # Verwende sigma-clipping für robuste Statistik
        _, median_val, _ = sigma_clipped_stats(annulus_data_1d, sigma=3.0)
        bkg_median.append(median_val)
    
    bkg_median = np.array(bkg_median)
    
    # Berechne die Gesamthintergrundkontamination in der Signal-Apertur
    # Fläche der Apertur = pi * r^2
    aperture_area = np.pi * aperture_radius**2
    total_bkg = bkg_median * aperture_area
    
    # Hintergrundkorrigierter Flux
    aperture_sum_bkgsub = phot_table['aperture_sum'] - total_bkg
    
    # Füge alle Informationen zur Tabelle hinzu
    phot_table['x'] = sources['xcentroid']
    phot_table['y'] = sources['ycentroid']
    phot_table['sky_median'] = bkg_median
    phot_table['aperture_sum_bkgsub'] = aperture_sum_bkgsub
    
    return phot_table, aperture_radius  # Gebe auch den verwendeten Radius zurück!


def detect_and_psf_photometry(data: np.ndarray, fwhm: float = 3.0, threshold_sigma: float = 5.0, aperture_radius: float = None):
    """Erkenne Quellen mit DAOStarFinder und führe PSF-Photometrie durch.
    
    PSF-Photometrie ist besser geeignet für dichte Sternfelder, da sie einzelne Sterne
    auch bei Überlappungen separieren kann.
    
    Liefert (phot_table, aperture_radius):
    - phot_table: Tabelle mit x, y, flux_fit (PSF-angepasster Flux)
    - aperture_radius: Der verwendete Aperturradius (für Kompatibilität mit Aperturphotometrie)
    
    Parameters:
        aperture_radius: Wenn None, wird automatisch 3.0×FWHM verwendet.
    
    Returns:
        (phot_table, aperture_radius): Tuple aus Photometrie-Tabelle und verwendetem Aperturradius
    """
    if not HAVE_PHOTUTILS:
        msg = (
            "Das Paket 'photutils' ist nicht installiert. Bitte installieren Sie es mit:\n"
            "pip install photutils\n"
            f"Original Fehler: {PHOTUTILS_ERROR}"
        )
        raise RuntimeError(msg)
    
    if not HAVE_PSF:
        print("  Warnung: PSF-Photometrie nicht verfügbar. Fallback auf Aperturphotometrie...")
        return detect_and_photometry(data, fwhm=fwhm, threshold_sigma=threshold_sigma, 
                                    aperture_radius=aperture_radius)
    
    # Automatische Apertur-Berechnung basierend auf FWHM
    if aperture_radius is None:
        aperture_radius = 3.0 * fwhm
        print(f"  Automatische Apertur: radius = 3.0 × FWHM = 3.0 × {fwhm:.1f} = {aperture_radius:.1f} Pixel")
    else:
        print(f"  Manuelle Apertur: radius = {aperture_radius:.1f} Pixel")
    
    # Robust stats für Quellenfindung
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    threshold = median + threshold_sigma * std
    
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)
    sources = daofind(data - median)
    if sources is None or len(sources) == 0:
        return None
    
    print(f"  PSF-Photometrie: Gefundene {len(sources)} Quellen, wende iterative PSF-Anpassung an...")
    
    try:
        from photutils.psf import IterativePSFPhotometry, GaussianPSF
        
        # Erstelle ein einfaches Gauss-PSF-Modell mit bekanntem FWHM
        psf_model = GaussianPSF(x_fwhm=fwhm, y_fwhm=fwhm)
        
        # Definiere die Größe der Anpassungsregion um jeden Stern
        fit_shape = (int(2*aperture_radius+1), int(2*aperture_radius+1))
        
        # Erstelle das PSF-Photometrie-Objekt
        photometry = IterativePSFPhotometry(
            psf_model=psf_model,
            fit_shape=fit_shape,
            finder=daofind,
            aperture_radius=aperture_radius,
            maxiters=3  # 3 Iterationen für iterative Subtraktion
        )
        
        # Führe PSF-Photometrie durch
        phot_table = photometry(data - median)
        
        if phot_table is None or len(phot_table) == 0:
            print("  Warnung: PSF-Photometrie lieferte keine Ergebnisse. Fallback auf Aperturphotometrie...")
            return detect_and_photometry(data, fwhm=fwhm, threshold_sigma=threshold_sigma, 
                                        aperture_radius=aperture_radius)
        
        # Stelle sicher, dass alle notwendigen Spalten existieren
        if 'x_fit' in phot_table.colnames:
            phot_table['x'] = phot_table['x_fit']
        else:
            phot_table['x'] = sources['xcentroid']
            
        if 'y_fit' in phot_table.colnames:
            phot_table['y'] = phot_table['y_fit']
        else:
            phot_table['y'] = sources['ycentroid']
        
        # Benenne Flux-Spalte um zur Kompatibilität
        if 'flux_fit' in phot_table.colnames:
            phot_table['aperture_sum_bkgsub'] = phot_table['flux_fit']
        elif 'flux' in phot_table.colnames:
            phot_table['aperture_sum_bkgsub'] = phot_table['flux']
        else:
            # Fallback: verwende die erste Flux-ähnliche Spalte
            for col in phot_table.colnames:
                if 'flux' in col.lower():
                    phot_table['aperture_sum_bkgsub'] = phot_table[col]
                    break
            else:
                # Wenn nichts gefunden, Aperturphotometrie als Fallback
                print("  Warnung: Keine Flux-Spalte in PSF-Ergebnissen. Fallback auf Aperturphotometrie...")
                return detect_and_photometry(data, fwhm=fwhm, threshold_sigma=threshold_sigma, 
                                            aperture_radius=aperture_radius)
        
        # Füge Himmelshintergrund hinzu (für Konsistenz)
        phot_table['sky_median'] = np.full(len(phot_table), median)
        
        print(f"  PSF-Photometrie erfolgreich: {len(phot_table)} Sterne angepasst")
        return phot_table, aperture_radius
        
    except Exception as e:
        print(f"  Warnung: PSF-Photometrie fehlgeschlagen ({type(e).__name__}: {e}). Fallback auf Aperturphotometrie...")
        # Fallback auf Aperturphotometrie
        positions = list(zip(sources['xcentroid'], sources['ycentroid']))
        apertures = CircularAperture(positions, r=aperture_radius)
        phot_table = aperture_photometry(data, apertures)
        
        # Hintergrund-Annulus
        r_in = aperture_radius * 1.5
        r_out = aperture_radius * 2.5
        annulus_apertures = CircularAnnulus(positions, r_in=r_in, r_out=r_out)
        annulus_masks = annulus_apertures.to_mask(method='center')
        
        bkg_median = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(data)
            annulus_data_1d = annulus_data[mask.data > 0]
            _, median_val, _ = sigma_clipped_stats(annulus_data_1d, sigma=3.0)
            bkg_median.append(median_val)
        
        bkg_median = np.array(bkg_median)
        aperture_area = np.pi * aperture_radius**2
        total_bkg = bkg_median * aperture_area
        aperture_sum_bkgsub = phot_table['aperture_sum'] - total_bkg
        
        phot_table['x'] = sources['xcentroid']
        phot_table['y'] = sources['ycentroid']
        phot_table['sky_median'] = bkg_median
        phot_table['aperture_sum_bkgsub'] = aperture_sum_bkgsub
        
        return phot_table, aperture_radius


def detect_and_photometry_wrapper(data: np.ndarray, method: str = 'aperture', fwhm: float = 3.0, 
                                   threshold_sigma: float = 5.0, aperture_radius: float = None):
    """Wrapper-Funktion, die zwischen Apertur- und PSF-Photometrie wählt.
    
    Parameters:
        data: 2D-Array mit Bilddaten
        method: 'aperture' oder 'psf'
        fwhm: FWHM in Pixel für Quellenfindung
        threshold_sigma: Schwelle in Sigma
        aperture_radius: Aperturradius (auto-berechnet falls None)
    
    Returns:
        (phot_table, aperture_radius): Photometrie-Ergebnisse
    """
    if method == 'psf':
        return detect_and_psf_photometry(data, fwhm=fwhm, threshold_sigma=threshold_sigma, 
                                         aperture_radius=aperture_radius)
    else:  # default: 'aperture'
        return detect_and_photometry(data, fwhm=fwhm, threshold_sigma=threshold_sigma, 
                                     aperture_radius=aperture_radius)


def estimate_transform(coords1, coords2, matches):
    """Schätze Translations-Offset zwischen zwei Koordinatensätzen basierend auf Matches."""
    if len(matches) == 0:
        return None
    
    offsets = []
    for i1, i2 in matches:
        offset = coords2[i2] - coords1[i1]
        offsets.append(offset)
    
    # Robuste Schätzung: Median-Offset
    median_offset = np.median(offsets, axis=0)
    return median_offset


def match_sources(tab1, tab2, max_sep: float = 5.0, search_radius: float = 100.0, 
                  filter_radius: Optional[float] = None, center: Optional[Tuple[float, float]] = None,
                  annulus_inner: Optional[float] = None):
    """Match sources in tab1 to tab2 using robust iterative pattern matching.
    
    Diese Methode ist robust gegenüber Verschiebungen zwischen den Bildern.
    Sie verwendet RANSAC-artiges iteratives Matching um die Transformation zu finden.
    
    WICHTIG: Die Radius-/Annulus-Filterung erfolgt NUR auf Basis von tab1 (B-Filter).
    Alle Matches werden dann nur für Sterne innerhalb dieses Bereichs in tab1 erstellt.
    
    Args:
        tab1, tab2: Photometry tables with 'x' and 'y' columns
        max_sep: Maximum separation for final matching (pixels)
        search_radius: Initial search radius for finding transformation (pixels)
        filter_radius: Optional outer radius - only sources within this radius from center in tab1
        center: Optional center point (x, y) for radius filter in tab1 (B-Filter) coordinates
        annulus_inner: Optional inner radius for annulus filter. When set, only sources between
                       annulus_inner and filter_radius are used.
    
    Rückgabe: (matched_pairs, offset) - Liste von Tupeln (index1, index2) und Offset [dx, dy]
    """
    if tab1 is None or tab2 is None:
        return [], None
    
    coords1 = np.vstack((np.array(tab1['x']), np.array(tab1['y']))).T
    coords2 = np.vstack((np.array(tab2['x']), np.array(tab2['y']))).T
    
    # Bestimme welche Sterne in tab1 (B-Filter) im Annulus oder Kreis liegen
    indices1_in_radius = None
    if filter_radius is not None and center is not None:
        dist1 = np.linalg.norm(coords1 - np.array(center), axis=1)
        if annulus_inner is not None:
            # Annulus-Filter: nur Sterne zwischen r_inner und r_outer (=filter_radius)
            mask1 = (dist1 >= annulus_inner) & (dist1 <= filter_radius)
            indices1_in_radius = set(np.where(mask1)[0])
            print(f"Annulus-Filter (B-Filter): {len(indices1_in_radius)}/{len(coords1)} Quellen im Ring "
                  f"({annulus_inner:.1f} – {filter_radius:.1f} Pixel) um ({center[0]:.1f}, {center[1]:.1f})")
        else:
            mask1 = dist1 <= filter_radius
            indices1_in_radius = set(np.where(mask1)[0])
            print(f"Radius-Filter (B-Filter): {len(indices1_in_radius)}/{len(coords1)} Quellen innerhalb von {filter_radius:.1f} Pixel vom Zentrum ({center[0]:.1f}, {center[1]:.1f})")
        
        if len(indices1_in_radius) < 3:
            print("Zu wenige Quellen im B-Filter nach Radius-/Annulus-Filterung (< 3).")
            return [], None
    
    if len(coords1) < 3 or len(coords2) < 3:
        print("Zu wenige Quellen für robustes Matching (< 3 in einem Bild).")
        return [], None
    
    # Sortiere nach Helligkeit für RANSAC (verwende alle Sterne für Offset-Bestimmung)
    if 'aperture_sum_bkgsub' in tab1.colnames:
        flux1 = np.array(tab1['aperture_sum_bkgsub'])
        flux2 = np.array(tab2['aperture_sum_bkgsub'])
        sort1 = np.argsort(flux1)[::-1]
        sort2 = np.argsort(flux2)[::-1]
        coords1_bright = coords1[sort1]
        coords2_bright = coords2[sort2]
        indices1_bright = sort1
        indices2_bright = sort2
    else:
        coords1_bright = coords1
        coords2_bright = coords2
        indices1_bright = np.arange(len(coords1))
        indices2_bright = np.arange(len(coords2))
    
    # Nimm die N hellsten Sterne für initiales Matching
    n_bright = min(100, len(coords1_bright), len(coords2_bright))
    
    best_offset = None
    best_n_matches = 0
    
    # RANSAC-artiger Ansatz: Teste verschiedene Hypothesen
    # Verwende bevorzugt Sterne innerhalb des Radius für die Offset-Bestimmung
    n_iterations = min(50, n_bright)
    
    for iteration in range(n_iterations):
        # Wähle einen Stern aus Bild 1, bevorzuge Sterne im Radius
        if iteration < min(10, n_bright):
            # Erste Iterationen: systematisch die hellsten
            i1_bright = iteration
        else:
            # Später: zufällige Auswahl, bevorzuge Sterne im Radius
            if indices1_in_radius is not None and len(indices1_in_radius) > 0 and np.random.random() < 0.7:
                # 70% Chance: wähle einen Stern im Radius
                candidates = [i for i in range(n_bright) if indices1_bright[i] in indices1_in_radius]
                if candidates:
                    i1_bright = np.random.choice(candidates)
                else:
                    i1_bright = np.random.randint(0, n_bright)
            else:
                i1_bright = np.random.randint(0, n_bright)
        
        star1 = coords1_bright[i1_bright]
        
        # Suche in Bild 2 nach möglichen Matches
        for i2_bright in range(min(n_bright, len(coords2_bright))):
            star2 = coords2_bright[i2_bright]
            
            # Hypothetischer Offset
            offset = star2 - star1
            
            # Teste diesen Offset: transformiere alle Sterne aus Bild 1
            coords1_transformed = coords1 + offset
            
            # Zähle wie viele Matches wir mit diesem Offset bekommen
            matches = []
            used2 = set()
            for j1 in range(len(coords1)):
                dists = np.linalg.norm(coords2 - coords1_transformed[j1], axis=1)
                j2 = int(np.argmin(dists))
                if dists[j2] <= max_sep and j2 not in used2:
                    matches.append((j1, j2))
                    used2.add(j2)
            
            # Ist dies die beste Hypothese bisher?
            if len(matches) > best_n_matches:
                best_n_matches = len(matches)
                best_offset = offset
                
                # Früher Abbruch wenn wir genug Matches haben
                if best_n_matches > min(len(coords1), len(coords2)) * 0.3:
                    break
        
        if best_n_matches > min(len(coords1), len(coords2)) * 0.3:
            break
    
    # Wende den besten Offset an und finde finale Matches
    if best_offset is not None:
        coords1_transformed = coords1 + best_offset
        
        # Finde alle Matches
        all_matches = []
        used2 = set()
        for i1 in range(len(coords1)):
            dists = np.linalg.norm(coords2 - coords1_transformed[i1], axis=1)
            i2 = int(np.argmin(dists))
            
            if dists[i2] <= max_sep and i2 not in used2:
                all_matches.append((i1, i2))
                used2.add(i2)
        
        # Filtere Matches: NUR Sterne die im B-Filter (tab1) innerhalb des Radius liegen
        if indices1_in_radius is not None:
            filtered_matches = [(i1, i2) for i1, i2 in all_matches if i1 in indices1_in_radius]
            print(f"Gefundener Offset zwischen Bildern: dx={best_offset[0]:.1f}, dy={best_offset[1]:.1f} Pixel")
            print(f"Matches gesamt: {len(all_matches)}, nach Radius-Filterung (B-Filter): {len(filtered_matches)}")
            return filtered_matches, best_offset
        else:
            print(f"Gefundener Offset zwischen Bildern: dx={best_offset[0]:.1f}, dy={best_offset[1]:.1f} Pixel")
            print(f"Matches nach Transformation: {len(all_matches)}")
            return all_matches, best_offset
    
    print("Warnung: Kein robuster Offset gefunden. Keine Matches möglich.")
    return [], None


def select_reference_star(images, phot_tables, matched_indices, aperture_radius, calibration_file=None, output_folder=None):
    """Interaktive Auswahl eines Referenzsterns und Eingabe seiner wahren Magnituden.
    
    Args:
        calibration_file: Optional - Pfad zum Speichern der Kalibrierungsdaten
        output_folder: Optional - Ordner für Ausgabedateien (PNG, etc.)
    
    Returns:
        (ref_index_B, mag_B_true, mag_V_true, ref_x, ref_y) oder (None, None, None, None, None) bei Abbruch
    """
    if len(images) < 1 or phot_tables[0] is None:
        print("Keine Daten für Referenzstern-Auswahl verfügbar.")
        return None, None, None, None, None
    
    data, header = images[0]  # B-Filter Bild
    tab = phot_tables[0]
    
    # Erstelle Plot für interaktive Auswahl
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_title('Klicken Sie auf einen Referenzstern im B-Filter Bild\n(Nur gematchte Sterne werden angezeigt)')
    
    try:
        norm = ImageNormalize(data, interval=ZScaleInterval())
    except Exception:
        norm = None
    
    ax.imshow(data, cmap='gray', origin='lower', norm=norm)
    
    # Speichere Koordinaten der gematchten Sterne für Click-Erkennung
    if matched_indices is not None and len(matched_indices) > 0:
        indices_B = [idx1 for idx1, idx2 in matched_indices]
        coords = np.vstack((np.array(tab['x']), np.array(tab['y']))).T
        coords_matched = coords[indices_B]
    else:
        coords_matched = np.array([])
        indices_B = []
    
    # Speichere Auswahl
    selected_ref = {'index': None, 'pos': None}
    
    def onclick(event):
        if event.inaxes != ax:
            return
        
        # Finde nächsten gematchten Stern
        click_pos = np.array([event.xdata, event.ydata])
        dists = np.linalg.norm(coords_matched - click_pos, axis=1)
        nearest_idx = np.argmin(dists)
        
        if dists[nearest_idx] < aperture_radius * 3:  # Innerhalb 3x Aperturradius
            selected_ref['index'] = indices_B[nearest_idx]
            selected_ref['pos'] = coords_matched[nearest_idx]
            
            # Markiere Auswahl
            ax.plot(selected_ref['pos'][0], selected_ref['pos'][1], 'r*', 
                   markersize=20, markeredgewidth=2, markeredgecolor='white')
            fig.canvas.draw()
            
            print(f"\nReferenzstern ausgewählt: Index {selected_ref['index']}, "
                  f"Position ({selected_ref['pos'][0]:.1f}, {selected_ref['pos'][1]:.1f})")
            print("Schließen Sie das Fenster um fortzufahren...")
    
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    fig.canvas.mpl_disconnect(cid)
    
    if selected_ref['index'] is None:
        print("Kein Referenzstern ausgewählt.")
        return None, None, None, None, None
    
    # Eingabe der wahren Magnituden
    print("\n=== Eingabe der wahren Magnituden des Referenzsterns ===")
    try:
        mag_B_true = float(input("B-Magnitude des Referenzsterns: "))
        mag_V_true = float(input("V-Magnitude des Referenzsterns: "))
    except (ValueError, EOFError):
        print("Ungültige Eingabe. Kalibrierung abgebrochen.")
        return None, None, None, None, None
    
    print(f"\nReferenzstern kalibriert: B={mag_B_true:.3f} mag, V={mag_V_true:.3f} mag")
    
    # Speichere Bild mit markiertem Referenzstern
    if output_folder:
        png_path = get_output_path(output_folder, 'reference_star.png')
        try:
            fig.savefig(png_path, dpi=150, bbox_inches='tight')
            print(f"Bild mit Referenzstern gespeichert: {png_path}")
        except Exception as e:
            print(f"WARNUNG: Konnte Bild nicht speichern: {e}", file=sys.stderr)
    
    # Speichere Kalibrierungsdaten falls Datei angegeben
    if calibration_file:
        save_calibration(calibration_file, selected_ref['index'], 
                        selected_ref['pos'][0], selected_ref['pos'][1],
                        mag_B_true, mag_V_true)
    
    return selected_ref['index'], mag_B_true, mag_V_true, selected_ref['pos'][0], selected_ref['pos'][1]


def save_calibration(filename, ref_index, ref_x, ref_y, mag_B, mag_V):
    """Speichere Kalibrierungsdaten in einer Textdatei."""
    try:
        with open(filename, 'w') as f:
            f.write("# Photometrie Kalibrierungsdatei\n")
            f.write(f"# Erstellt am: {np.datetime64('now')}\n")
            f.write("#\n")
            f.write("# Format: ref_index ref_x ref_y mag_B mag_V\n")
            f.write(f"{ref_index} {ref_x:.6f} {ref_y:.6f} {mag_B:.6f} {mag_V:.6f}\n")
        print(f"Kalibrierungsdaten gespeichert in: {filename}")
    except Exception as e:
        print(f"WARNUNG: Konnte Kalibrierungsdaten nicht speichern: {e}", file=sys.stderr)


def save_photometry_table(filename, mags_V, colors_BV, x_coords=None, y_coords=None):
    """Speichere Photometrie-Tabelle mit V und B-V Werten.
    
    Args:
        filename: Pfad zur Ausgabedatei
        mags_V: Array mit V-Magnituden
        colors_BV: Array mit B-V Farben
        x_coords: Optional - x-Koordinaten der Sterne
        y_coords: Optional - y-Koordinaten der Sterne
    """
    try:
        with open(filename, 'w') as f:
            f.write("# Photometrie-Tabelle des Sternhaufens\n")
            f.write(f"# Erstellt am: {np.datetime64('now')}\n")
            f.write("#\n")
            if x_coords is not None and y_coords is not None:
                f.write("# Format: index x y V B-V\n")
                f.write("# x y: Pixel-Koordinaten im B-Filter Bild\n")
            else:
                f.write("# Format: index V B-V\n")
            f.write("# V: V-Magnitude (mag)\n")
            f.write("# B-V: Farbindex (mag)\n")
            f.write("#\n")
            
            for i in range(len(mags_V)):
                if x_coords is not None and y_coords is not None:
                    f.write(f"{i:4d} {x_coords[i]:8.2f} {y_coords[i]:8.2f} {mags_V[i]:7.3f} {colors_BV[i]:7.3f}\n")
                else:
                    f.write(f"{i:4d} {mags_V[i]:7.3f} {colors_BV[i]:7.3f}\n")
        
        print(f"Photometrie-Tabelle gespeichert: {filename}")
        print(f"  {len(mags_V)} Sterne gespeichert")
    except Exception as e:
        print(f"WARNUNG: Konnte Photometrie-Tabelle nicht speichern: {e}", file=sys.stderr)


def load_calibration(filename, phot_tables, matched_indices):
    """Lade Kalibrierungsdaten aus einer Textdatei.
    
    Returns:
        (ref_index_B, mag_B_true, mag_V_true, ref_x, ref_y) oder (None, None, None, None, None) bei Fehler
    """
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                
                parts = line.split()
                if len(parts) >= 5:
                    ref_index = int(parts[0])
                    ref_x = float(parts[1])
                    ref_y = float(parts[2])
                    mag_B = float(parts[3])
                    mag_V = float(parts[4])
                    
                    print(f"Kalibrierungsdaten geladen aus: {filename}")
                    print(f"  Referenzstern: Index {ref_index}, Position ({ref_x:.1f}, {ref_y:.1f})")
                    print(f"  Magnituden: B={mag_B:.3f}, V={mag_V:.3f}")
                    
                    # Prüfe ob der Index in den gematchten Sternen existiert
                    if phot_tables and phot_tables[0] is not None:
                        if ref_index >= len(phot_tables[0]):
                            print(f"WARNUNG: Referenzstern-Index {ref_index} existiert nicht in aktuellen Daten!")
                            return None, None, None, None, None
                    
                    # Prüfe ob der Stern in matched_indices ist
                    if matched_indices:
                        found = any(idx1 == ref_index for idx1, idx2 in matched_indices)
                        if not found:
                            print(f"WARNUNG: Referenzstern ist nicht in den gematchten Sternen!")
                            print(f"  Möglicherweise hat sich die Sternauswahl geändert (anderer Radius?).")
                            return None, None, None, None, None
                    
                    return ref_index, mag_B, mag_V, ref_x, ref_y
        
        print(f"FEHLER: Keine gültigen Kalibrierungsdaten in {filename} gefunden.")
        return None, None, None, None, None
    
    except FileNotFoundError:
        print(f"FEHLER: Kalibrierungsdatei nicht gefunden: {filename}")
        return None, None, None, None, None
    except Exception as e:
        print(f"FEHLER beim Laden der Kalibrierungsdatei: {e}", file=sys.stderr)
        return None, None, None, None, None


def calibrate_magnitudes(matched, tabB, tabV, ref_index_B, mag_B_true, mag_V_true):
    """Kalibriere alle Magnituden basierend auf einem Referenzstern.
    
    Returns:
        (calibrated_mags_B, calibrated_mags_V, colors_BV, matched_indices_valid)
    """
    sumsB = np.asarray(tabB['aperture_sum_bkgsub'])
    sumsV = np.asarray(tabV['aperture_sum_bkgsub'])
    
    # Finde den Referenzstern in den Matches
    ref_match_idx = None
    for i, (iB, iV) in enumerate(matched):
        if iB == ref_index_B:
            ref_match_idx = i
            ref_index_V = iV
            break
    
    if ref_match_idx is None:
        print("FEHLER: Referenzstern nicht in den Matches gefunden!")
        return None, None, None, None
    
    # Berechne instrumentale Magnituden des Referenzsterns
    ref_flux_B = sumsB[ref_index_B]
    ref_flux_V = sumsV[ref_index_V]
    
    if ref_flux_B <= 0 or ref_flux_V <= 0:
        print("FEHLER: Referenzstern hat negativen oder Null-Flux!")
        return None, None, None, None
    
    mag_B_inst_ref = -2.5 * np.log10(ref_flux_B)
    mag_V_inst_ref = -2.5 * np.log10(ref_flux_V)
    
    # Berechne Nullpunkte
    zeropoint_B = mag_B_true - mag_B_inst_ref
    zeropoint_V = mag_V_true - mag_V_inst_ref
    
    print(f"\n=== Kalibrierung Debug ===")
    print(f"Referenzstern: B_true = {mag_B_true:.3f}, V_true = {mag_V_true:.3f}, (B-V)_true = {mag_B_true - mag_V_true:.3f}")
    print(f"Instrumentale Magnituden: B_inst = {mag_B_inst_ref:.3f}, V_inst = {mag_V_inst_ref:.3f}")
    print(f"Nullpunkte: ZP_B = {zeropoint_B:.3f} mag, ZP_V = {zeropoint_V:.3f} mag")
    print(f"Differenz: ZP_V - ZP_B = {zeropoint_V - zeropoint_B:.3f} mag")
    
    # Kalibriere alle Sterne
    mags_B = []
    mags_V = []
    colors = []
    valid_matches = []
    
    for iB, iV in matched:
        flux_B = sumsB[iB]
        flux_V = sumsV[iV]
        
        if flux_B > 0 and flux_V > 0:
            mag_B = -2.5 * np.log10(flux_B) + zeropoint_B
            mag_V = -2.5 * np.log10(flux_V) + zeropoint_V
            color_BV = mag_B - mag_V
            
            mags_B.append(mag_B)
            mags_V.append(mag_V)
            colors.append(color_BV)
            valid_matches.append((iB, iV))
    
    colors_array = np.array(colors)
    print(f"\nKalibrierte Photometrie:")
    print(f"  {len(colors_array)} Sterne erfolgreich kalibriert")
    print(f"  B-V: min={np.min(colors_array):.3f}, max={np.max(colors_array):.3f}, median={np.median(colors_array):.3f}")
    print(f"  WARNUNG: Wenn B-V median < 0, könnte B und V vertauscht sein!")
    
    return np.array(mags_B), np.array(mags_V), colors_array, valid_matches


def instrumental_magnitude(flux: float, zero_point: float = 25.0) -> Optional[float]:
    """Berechne instrumentale Magnitude aus Flux.
    
    Formel: mag = -2.5 * log10(flux) + zero_point
    
    Der zero_point sollte so gewählt werden, dass typische Sterne positive Magnituden haben.
    Standard: 25.0 (ergibt positive Magnituden für typische Flüsse)
    """
    try:
        if flux <= 0:
            return None
        return -2.5 * np.log10(flux) + zero_point
    except Exception:
        return None


def load_parsec_table(filepath, age_log, main_sequence_only=False):
    """Lade PARSEC-Isochrone aus einer Multi-Age-Tabelle.
    
    Diese Funktion liest CMD 3.7/3.8 Output-Dateien die mehrere Alter enthalten.
    
    Args:
        filepath: Pfad zur PARSEC-Tabelle (z.B. isochrones.txt)
        age_log: Gewünschtes log(age) in Jahren
        main_sequence_only: Wenn True, nur Hauptreihen-Sterne (label <= 1)
    
    Returns:
        (B-V, M_V) Arrays oder (None, None) bei Fehler
    """
    print(f"  Lade PARSEC-Tabelle: {filepath}")
    
    try:
        # Lese alle Zeilen
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Finde Header-Zeile um Spalten zu identifizieren
        header_line = None
        data_start = 0
        for i, line in enumerate(lines):
            if line.startswith('# Zini'):
                header_line = line
                data_start = i + 1
                break
        
        if header_line is None:
            print(f"  WARNUNG: Kein Header gefunden in {filepath}")
            return None, None
        
        # Parse Header um Spalten-Indices zu finden
        header_parts = header_line.replace('#', '').strip().split()
        
        try:
            col_logage = header_parts.index('logAge')
            col_bmag = header_parts.index('Bmag')
            col_vmag = header_parts.index('Vmag')
            col_label = header_parts.index('label') if 'label' in header_parts else None
            if main_sequence_only and col_label is not None:
                print(f"  Spalten gefunden: logAge={col_logage}, Bmag={col_bmag}, Vmag={col_vmag}, label={col_label} (nur Hauptreihe)")
            else:
                print(f"  Spalten gefunden: logAge={col_logage}, Bmag={col_bmag}, Vmag={col_vmag}")
        except ValueError as e:
            print(f"  FEHLER: Spalte nicht gefunden: {e}")
            return None, None
        
        # Sammle Daten für das gewünschte Alter
        b_v_list = []
        v_abs_list = []
        ages_found = set()
        n_total = 0
        n_filtered = 0
        
        tolerance = 0.05  # Toleranz für log(age) Vergleich
        
        for line in lines[data_start:]:
            # Überspringe Kommentare und leere Zeilen
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.split()
            if len(parts) < max(col_logage, col_bmag, col_vmag) + 1:
                continue
            
            try:
                logage = float(parts[col_logage])
                ages_found.add(logage)
                
                # Prüfe ob dieses Alter passt
                if abs(logage - age_log) < tolerance:
                    n_total += 1
                    
                    # Filtere nach evolutionärem Stadium wenn gewünscht
                    if main_sequence_only and col_label is not None:
                        try:
                            label = int(parts[col_label])
                            # PARSEC labels: 0 = pre-MS, 1 = MS, 2 = SGB, 3 = RGB, etc.
                            # Wir wollen nur 0 und 1 (pre-MS und MS)
                            if label > 1:
                                n_filtered += 1
                                continue
                        except (ValueError, IndexError):
                            pass
                    
                    b_mag = float(parts[col_bmag])
                    v_mag = float(parts[col_vmag])
                    
                    # Filtere unphysikalische Werte
                    if -15 < v_mag < 20 and -1 < (b_mag - v_mag) < 3:
                        b_v_list.append(b_mag - v_mag)
                        v_abs_list.append(v_mag)
            except (ValueError, IndexError):
                continue
        
        if len(b_v_list) == 0:
            print(f"  WARNUNG: Keine Daten für log(age)={age_log:.3f} gefunden!")
            if ages_found:
                available = sorted(list(ages_found))[:10]  # Zeige erste 10
                print(f"  Verfügbare Alter: {', '.join([f'{a:.2f}' for a in available])}...")
            return None, None
        
        b_v = np.array(b_v_list)
        v_abs = np.array(v_abs_list)
        
        if main_sequence_only and n_filtered > 0:
            print(f"  Erfolgreich! {len(b_v)} Hauptreihen-Sterne (gefiltert: {n_filtered} Post-MS)")
        else:
            print(f"  Erfolgreich! {len(b_v)} Datenpunkte für log(age)={age_log:.3f}")
        print(f"  B-V: [{np.min(b_v):.2f}, {np.max(b_v):.2f}], M_V: [{np.min(v_abs):.2f}, {np.max(v_abs):.2f}]")
        
        return b_v, v_abs
        
    except Exception as e:
        print(f"  Fehler beim Laden: {e}")
        return None, None


def download_parsec_isochrone_cmd(age_log, z=0.0152):
    """Lade echte PARSEC-Isochrone vom CMD 3.7 Web-Service.
    
    HINWEIS: Automatischer Download kann instabil sein!
    Für beste Ergebnisse:
    1. Gehen Sie zu: http://stev.oapd.inaf.it/cgi-bin/cmd
    2. Wählen Sie:
       - Evolutionary tracks: PARSEC v1.2S
       - Photometric system: Bessell (Johnson-Cousins)
       - Output: Isochrone tables
       - log(age/yr): Ihr gewünschtes Alter (z.B. 7.0 für 10 Myr)
       - Metallicity Z: 0.0152 (solar)
    3. Klicken Sie "Submit"
    4. Laden Sie die Datei herunter
    5. Speichern Sie als 'isochrone_<log_age>.dat' im FITS-Ordner
       ODER extrahieren Sie Spalten Bmag und Vmag und speichern als 2-Spalten-Datei (B-V, M_V)
    
    Args:
        age_log: Logarithmus des Alters in Jahren (z.B. 7.0 für 10 Millionen Jahre)
        z: Metallizität (Standard: 0.0152 = solar)
    
    Returns:
        (B-V, M_V) Arrays oder (None, None) bei Fehler
    """
    import urllib.request
    import urllib.parse
    import time
    
    # CMD 3.7 Web-Service URL
    base_url = "http://stev.oapd.inaf.it/cgi-bin/cmd"
    
    # Parameter für PARSEC v1.2S Isochrone
    params = {
        'cmd_version': '3.7',
        'track_parsec': 'parsec_CAF09_v1.2S',
        'track_colibri': 'parsec_CAF09_v1.2S_S_LMC_08_web',
        'track_postagb': 'no',
        'n_inTPC': '10',
        'eta_reimers': '0.2',
        'kind_interp': '1',
        'kind_postagb': '-1',
        'photsys_file': 'tab_mag_odfnew/tab_mag_bessell.dat',  # Johnson-Cousins UBVRI
        'photsys_version': 'yang',
        'dust_sourceM': 'dpmod60alox40',
        'dust_sourceC': 'AMCSIC15',
        'kind_mag': '2',  # Absolute magnitudes
        'kind_dust': '0',  # No dust
        'extinction_av': '0.0',
        'imf_file': 'tab_imf/imf_kroupa_orig.dat',
        'isoc_isagelog': '1',  # log(age)
        'isoc_lagelow': str(age_log),
        'isoc_lageupp': str(age_log),
        'isoc_dlage': '0.0',
        'isoc_ismetlog': '0',  # Linear Z
        'isoc_zlow': str(z),
        'isoc_zupp': str(z),
        'isoc_dz': '0.0',
        'output_kind': '0',  # Isochrone table
        'output_evstage': '1',
        'lf_maginf': '-15',
        'lf_magsup': '20',
        'lf_deltamag': '0.5',
        'sim_mtot': '1.0e4',
        'submit_form': 'Submit'
    }
    
    try:
        print(f"  Fordere PARSEC-Isochrone an (log(age)={age_log:.2f}, Z={z:.4f})...")
        
        # POST-Request
        data = urllib.parse.urlencode(params).encode('ascii')
        req = urllib.request.Request(base_url, data=data)
        
        with urllib.request.urlopen(req, timeout=30) as response:
            result = response.read().decode('utf-8')
        
        # Extrahiere Download-URL aus der Antwort
        # Der Server gibt eine HTML-Seite mit einem Link zur Datei zurück
        import re
        match = re.search(r'(output\d+\.dat)', result)
        
        if not match:
            print("  WARNUNG: Keine Download-URL in Antwort gefunden")
            return None, None
        
        output_file = match.group(1)
        download_url = f"http://stev.oapd.inaf.it/~webmaster/tmp/{output_file}"
        
        print(f"  Lade Daten von: {download_url}")
        time.sleep(1)  # Kurze Pause um Server nicht zu überlasten
        
        with urllib.request.urlopen(download_url, timeout=30) as response:
            data_text = response.read().decode('utf-8')
        
        # Parse die Daten
        lines = data_text.strip().split('\n')
        b_v_list = []
        v_abs_list = []
        
        for line in lines:
            # Überspringe Kommentare und Header
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.split()
            if len(parts) < 20:
                continue
            
            try:
                # PARSEC Spalten (typisch):
                # Cols: Zini logAge Mini Mass logL logTe [B-V] U B V R I J H K
                # Die genaue Position hängt vom Output ab
                # Für bessell.dat: Spalten sind etwa bei Index 14-19
                
                # Typischerweise:
                # col 0: Z, 1: log(age), 2: M_ini, 3: Mass, 4: logL, 5: logTe
                # col 6-10: stage info
                # col ~14: B mag, ~15: V mag
                
                # Sicherer: Suche nach B und V in den letzten Spalten
                # Format: ... Bmag Vmag Rmag Imag ...
                if len(parts) >= 16:
                    b_mag = float(parts[14])  # B (absolute)
                    v_mag = float(parts[15])  # V (absolute)
                    
                    # Filtere unphysikalische Werte
                    if -15 < v_mag < 20 and -1 < (b_mag - v_mag) < 3:
                        b_v_list.append(b_mag - v_mag)
                        v_abs_list.append(v_mag)
            except (ValueError, IndexError):
                continue
        
        if len(b_v_list) == 0:
            print("  WARNUNG: Keine gültigen Datenpunkte gefunden")
            return None, None
        
        b_v = np.array(b_v_list)
        v_abs = np.array(v_abs_list)
        
        print(f"  Erfolgreich! {len(b_v)} Datenpunkte, B-V: [{np.min(b_v):.2f}, {np.max(b_v):.2f}]")
        
        return b_v, v_abs
        
    except Exception as e:
        print(f"  Fehler beim Laden von CMD 3.7: {e}")
        return None, None


def download_isochrone_cmd(age_log, z=0.0152, output_file=None, main_sequence_only=False):
    """Lade Isochrone - versuche verschiedene Methoden.
    
    Args:
        age_log: Logarithmus des Alters in Jahren (z.B. 7.0 für 10 Millionen Jahre)
        z: Metallizität (Standard: 0.0152 = solar)
        output_file: Optional - Pfad zum Speichern der Isochrone
        main_sequence_only: Wenn True, nur Hauptreihen-Sterne laden
    
    Returns:
        (B-V, V) Arrays oder (None, None) bei Fehler
    """
    print(f"Lade Isochrone für log(age)={age_log:.2f}, Z={z:.4f}...")
    
    # Methode 1: Suche nach isochrones.txt (Multi-Age PARSEC-Tabelle) im gleichen Ordner wie photo.py
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parsec_table = os.path.join(script_dir, 'isochrones.txt')
    
    if os.path.exists(parsec_table):
        print(f"  Gefunden: Multi-Age PARSEC-Tabelle")
        b_v, v_abs = load_parsec_table(parsec_table, age_log, main_sequence_only=main_sequence_only)
        if b_v is not None:
            return b_v, v_abs
        print(f"  Konnte log(age)={age_log:.2f} nicht in Tabelle finden, versuche andere Methoden...")
    
    # Methode 2: Versuche lokale Datei zu laden falls vorhanden (gespeicherte Einzelalter-Isochrone)
    if output_file and os.path.exists(output_file):
        print(f"  Lade gespeicherte Isochrone: {output_file}")
        try:
            # Versuche verschiedene Formate
            # Format 1: Einfach (B-V, M_V)
            data = np.loadtxt(output_file, comments='#')
            if len(data) > 0 and data.shape[1] >= 2:
                print(f"  Erfolgreich! {len(data)} Datenpunkte geladen")
                return data[:, 0], data[:, 1]
        except Exception as e:
            print(f"  Fehler beim Laden als 2-Spalten-Datei: {e}")
            
            # Format 2: PARSEC-Output (viele Spalten)
            try:
                b_v_list = []
                v_abs_list = []
                
                with open(output_file, 'r') as f:
                    for line in f:
                        if line.startswith('#') or not line.strip():
                            continue
                        
                        parts = line.split()
                        if len(parts) < 2:
                            continue
                        
                        try:
                            # Versuche verschiedene Spalten-Konfigurationen
                            if len(parts) >= 16:
                                # Voller PARSEC-Output
                                b_mag = float(parts[14])
                                v_mag = float(parts[15])
                            elif len(parts) >= 2:
                                # Vereinfachtes Format: nur B-V und M_V
                                b_v_val = float(parts[0])
                                v_mag = float(parts[1])
                                b_v_list.append(b_v_val)
                                v_abs_list.append(v_mag)
                                continue
                            else:
                                continue
                            
                            # Für volle PARSEC-Dateien
                            if -15 < v_mag < 20 and -1 < (b_mag - v_mag) < 3:
                                b_v_list.append(b_mag - v_mag)
                                v_abs_list.append(v_mag)
                        except (ValueError, IndexError):
                            continue
                
                if len(b_v_list) > 0:
                    print(f"  Erfolgreich! {len(b_v_list)} Datenpunkte aus PARSEC-Datei")
                    return np.array(b_v_list), np.array(v_abs_list)
            except Exception as e2:
                print(f"  Fehler beim Parsen: {e2}")
    
    # NICHT mehr versuchen vom Server zu laden - zu instabil!
    # Versuche echte PARSEC-Isochrone vom CMD 3.7 zu laden
    # b_v, v_mag = download_parsec_isochrone_cmd(age_log, z)
    # 
    # # Falls erfolgreich, speichere für zukünftige Verwendung
    # if b_v is not None and output_file:
    #     try:
    #         np.savetxt(output_file, np.column_stack([b_v, v_mag]), 
    #                   header=f'B-V V_abs (log(age)={age_log:.2f}, Z={z:.4f})')
    #         print(f"  Isochrone gespeichert: {output_file}")
    #     except Exception as e:
    #         print(f"  Warnung: Konnte nicht speichern: {e}")
    #     return b_v, v_mag
    
    # Fallback: Verwende theoretisches Modell
    print(f"  Konnte log(age)={age_log:.2f} nicht in isochrones.txt finden.")
    print("  Fallback: Verwende theoretisches Modell...")
    b_v, v_mag = generate_parsec_like_isochrone(age_log, z)
    
    # Speichere auch das theoretische Modell
    if output_file and b_v is not None:
        try:
            np.savetxt(output_file, np.column_stack([b_v, v_mag]), 
                      header=f'B-V V_abs (theoretical, log(age)={age_log:.2f}, Z={z:.4f})')
            print(f"  Theoretische Isochrone gespeichert: {output_file}")
        except Exception as e:
            print(f"  Warnung: Konnte nicht speichern: {e}")
    
    return b_v, v_mag


def generate_parsec_like_isochrone(age_log, z=0.0152):
    """Erzeuge realistische Isochrone basierend auf PARSEC-Modellen.
    
    Verwendet empirische Fits an PARSEC-Isochrone für verschiedene Alter.
    
    Args:
        age_log: Logarithmus des Alters in Jahren
        z: Metallizität (Standard: 0.0152 = solar)
    
    Returns:
        (B-V, V) Arrays
    """
    age = 10**age_log  # Jahre
    age_gyr = age / 1e9  # Gigajahre
    
    # PARSEC Haupt-Sequenz Relation (empirisch angepasst)
    # Basiert auf solarer Metallizität PARSEC v1.2S Isochrone
    
    # Massenbereiche für Hauptreihe
    # Wir modellieren Sterne von ~0.5 bis ~50 Sonnenmassen (relevant für Haufen)
    # MEHR Punkte für glattere Isochrone!
    masses = np.concatenate([
        np.linspace(0.5, 2.0, 100),    # Niedrige bis mittlere Massen
        np.linspace(2.0, 10.0, 200),   # Mittlere Massen (sehr gut aufgelöst)
        np.linspace(10.0, 50.0, 200)   # Hohe Massen (wichtig für junge Haufen!)
    ])
    
    # Turn-off Masse hängt vom Alter ab
    # Haupt-Sequenz Lifetime: t_MS ≈ 10 Gyr * (M/M_sun)^(-2.5)
    # NGC 869: ~13 Myr → M_turnoff ~ 8-10 M_sun
    
    if age_gyr < 0.001:  # < 1 Million Jahre - extrem jung
        m_turnoff = 50.0  # Alle Massen auf Hauptreihe
    else:
        # M_turnoff aus t_MS = age
        m_turnoff = (10.0 / age_gyr) ** 0.4
        m_turnoff = np.clip(m_turnoff, 0.8, 50.0)
    
    print(f"  Alter: {age/1e6:.1f} Myr, Turn-off Masse: {m_turnoff:.2f} M_sun")
    
    b_v_list = []
    v_abs_list = []
    
    for m in masses:
        if m <= m_turnoff:
            # Hauptreihe - verbesserte Farb-Temperatur-Relation
            # Basierend auf Spektraltypen und PARSEC-Daten
            # Glatte Interpolation für bessere Isochrone
            if m < 0.7:
                bv = 1.4  # M-Sterne (sehr rot)
            elif m < 1.0:
                # K-Sterne - sanfter Übergang
                t = (m - 0.7) / 0.3
                bv = 1.4 * (1-t) + 0.8 * t
            elif m < 1.5:
                # G-Sterne
                t = (m - 1.0) / 0.5
                bv = 0.8 * (1-t) + 0.65 * t
            elif m < 2.5:
                # F-Sterne
                t = (m - 1.5) / 1.0
                bv = 0.65 * (1-t) + 0.3 * t
            elif m < 7.0:
                # A-Sterne
                t = (m - 2.5) / 4.5
                bv = 0.3 * (1-t) + 0.0 * t
            elif m < 17.0:
                # B-Sterne (wichtig!)
                t = (m - 7.0) / 10.0
                bv = 0.0 * (1-t) + (-0.15) * t
            else:
                # O-Sterne
                t = min((m - 17.0) / 33.0, 1.0)
                bv = -0.15 * (1-t) + (-0.25) * t
            
            bv = np.clip(bv, -0.33, 1.5)
            
            # Absolute V-Magnitude - verbesserte Masse-Leuchtkraft-Relation
            # Wichtig: Für massive Sterne (M > 2 M_sun) ist die Relation flacher!
            if m < 0.43:
                # Sehr niedrige Massen: flache Relation
                log_l = 2.3 * np.log10(m)
            elif m < 2.0:
                # Sonnenähnliche Sterne: L ~ M^4
                log_l = 4.0 * np.log10(m)
            elif m < 20.0:
                # Massive Sterne: L ~ M^3.0 (flacher!)
                # Wichtig für korrekte absolute Helligkeit
                log_l = 3.0 * np.log10(m) + 0.3  # Offset für bessere Anpassung
            else:
                # Sehr massive Sterne: L ~ M^2.5
                log_l = 2.5 * np.log10(m) + 0.6
            
            v_abs = 4.83 - 2.5 * log_l
            
        else:
            # Post-Hauptreihe (Riesen, falls Turn-off erreicht)
            # Werden röter und heller
            bv = 0.8 + 0.3 * (m - m_turnoff) / m_turnoff
            bv = np.clip(bv, 0.5, 2.0)
            
            # Heller als Turn-off-Punkt
            m_to_log_l = 3.0 * np.log10(m_turnoff) + 0.3
            m_to_v = 4.83 - 2.5 * m_to_log_l
            v_abs = m_to_v - 2.0 * (m - m_turnoff) / m_turnoff
            v_abs = np.clip(v_abs, -8.0, m_to_v)
        
        b_v_list.append(bv)
        v_abs_list.append(v_abs)
    
    # Sortiere nach B-V für schöne Darstellung
    b_v_array = np.array(b_v_list)
    v_abs_array = np.array(v_abs_list)
    
    sort_idx = np.argsort(b_v_array)
    b_v_array = b_v_array[sort_idx]
    v_abs_array = v_abs_array[sort_idx]
    
    # Entferne Duplikate in B-V
    unique_bv = []
    unique_v = []
    for i in range(len(b_v_array)):
        if i == 0 or abs(b_v_array[i] - unique_bv[-1]) > 0.001:
            unique_bv.append(b_v_array[i])
            unique_v.append(v_abs_array[i])
    
    print(f"  Isochrone erstellt: {len(unique_bv)} Punkte, B-V: [{min(unique_bv):.2f}, {max(unique_bv):.2f}]")
    
    return np.array(unique_bv), np.array(unique_v)


def fit_isochrone_distance_age(obs_bv, obs_v, ebv, age_range=(6.5, 8.5), age_steps=20, output_folder=None, fixed_age_myr=None):
    """Bestimme Entfernung und Alter durch Isochrone-Anpassung.
    
    Args:
        obs_bv: Beobachtete B-V Farben
        obs_v: Beobachtete V-Magnituden
        ebv: Rötung E(B-V)
        age_range: (min_log_age, max_log_age) in Jahren (ignoriert wenn fixed_age_myr gesetzt)
        age_steps: Anzahl der zu testenden Alter (ignoriert wenn fixed_age_myr gesetzt)
        output_folder: Optional - Ordner zum Speichern von Plots
        fixed_age_myr: Optional - Fixiertes Alter in Millionen Jahren (z.B. 10.6)
    
    Returns:
        Dictionary mit: best_age_log, dm, distance_pc, chi2, iso_bv, iso_v, ebv, obs_bv_dereddened, obs_v_dereddened
        oder None bei Fehler
    """
    print(f"\n=== Isochrone-Anpassung ===")
    print(f"Rötung E(B-V) = {ebv:.3f} mag")
    
    if fixed_age_myr is not None:
        # Fixiertes Alter
        fixed_age_log = np.log10(fixed_age_myr * 1e6)
        print(f"Fixiertes Alter: {fixed_age_myr} Myr (log(age) = {fixed_age_log:.3f})")
        ages_log = np.array([fixed_age_log])
        age_steps = 1
    else:
        print(f"Altersbereich: log(age) = {age_range[0]:.2f} - {age_range[1]:.2f}")
        ages_log = np.linspace(age_range[0], age_range[1], age_steps)
    
    # Korrigiere Beobachtungen für Rötung
    # A_V = 3.1 * E(B-V) (Standard-Extinktionsgesetz nach Cardelli et al.)
    # A_B = 4.1 * E(B-V)
    av = 3.1 * ebv
    ab = 4.1 * ebv
    
    # Dereddening:
    # (B-V)_0 = (B-V)_obs - E(B-V)
    # V_0 = V_obs - A_V
    obs_bv_dereddened = obs_bv - ebv
    obs_v_dereddened = obs_v - av
    
    print(f"Extinktion A_V = {av:.3f} mag, A_B = {ab:.3f} mag")
    print(f"Beobachtete Sterne:")
    print(f"  B-V (beobachtet): min={np.min(obs_bv):.2f}, max={np.max(obs_bv):.2f}, median={np.median(obs_bv):.2f}")
    print(f"  (B-V)_0 (dereddened): min={np.min(obs_bv_dereddened):.2f}, max={np.max(obs_bv_dereddened):.2f}, median={np.median(obs_bv_dereddened):.2f}")
    print(f"  V (beobachtet): min={np.min(obs_v):.2f}, max={np.max(obs_v):.2f}, median={np.median(obs_v):.2f}")
    print(f"  V_0 (dereddened): min={np.min(obs_v_dereddened):.2f}, max={np.max(obs_v_dereddened):.2f}, median={np.median(obs_v_dereddened):.2f}")
    
    # Teste verschiedene Alter (oder nur ein fixiertes Alter)
    chi2_values = []
    distance_moduli = []
    
    if fixed_age_myr is not None:
        print(f"\nBestimme Entfernung für fixiertes Alter...")
    else:
        print(f"\nTeste {age_steps} Alter...")
    
    best_overall_age = None
    best_overall_dm = None
    best_overall_chi2 = np.inf
    
    for idx_age, age_log in enumerate(ages_log):
        # Lade Isochrone - NUR HAUPTREIHE für Chi²-Anpassung!
        iso_bv, iso_v = download_isochrone_cmd(age_log, main_sequence_only=True)
        
        if iso_bv is None or len(iso_bv) == 0:
            chi2_values.append(np.inf)
            distance_moduli.append(np.nan)
            continue
        
        # WICHTIG: Beschränke Isochrone auf den B-V Bereich der Daten!
        # Verwende BEOBACHTETE B-V (nicht dereddened) für Trimming!
        # Die Isochrone ist intrinsisch, wird dann gerötet zum Vergleich
        bv_min_data = np.min(obs_bv)  # BEOBACHTET
        bv_max_data = np.max(obs_bv)  # BEOBACHTET
        
        # Röte die Isochrone für den Vergleich
        iso_bv_reddened = iso_bv + ebv
        
        # Verwende nur Isochrone im BEOBACHTETEN Datenbereich (mit kleinem Puffer)
        iso_mask = (iso_bv_reddened >= bv_min_data - 0.1) & (iso_bv_reddened <= bv_max_data + 0.1)
        iso_bv_trimmed = iso_bv[iso_mask]  # Intrinsisch (für Fitting)
        iso_v_trimmed = iso_v[iso_mask]
        
        if len(iso_bv_trimmed) == 0:
            chi2_values.append(np.inf)
            distance_moduli.append(np.nan)
            continue
        
        # Finde bestes Entfernungsmodul durch vertikale Verschiebung
        # Minimiere Chi² über alle möglichen Verschiebungen
        # Für 1300 pc: DM = 5*log10(1300) - 5 = 10.57 mag
        # Teste größeren Bereich: 9 bis 14 mag (500 pc bis 6300 pc)
        
        dm_range = np.arange(9.0, 14.0, 0.05)  # Feinere Schritte für bessere Genauigkeit
        best_chi2 = np.inf
        best_dm = np.nan
        
        # Speichere Chi² Kurve für Debug (nur bei fixiertem Alter)
        chi2_curve = [] if fixed_age_myr is not None else None
        dm_curve = [] if fixed_age_myr is not None else None
        
        for dm in dm_range:
            # Verschiebe Isochrone (absolute V-Mag -> scheinbare V-Mag)
            # Verwende nur den getrimmten Teil der Isochrone!
            iso_v_shifted = iso_v_trimmed + dm
            
            # NEUE STRATEGIE: Anpassung an den BLAUEN RAND (ZAMS)
            # Für jeden Stern: Finde hellsten Isochrone-Punkt bei ähnlicher Farbe
            # Die Hauptreihe ist der blaue/helle Rand der Verteilung!
            
            chi2 = 0.0
            n_matched = 0
            total_weight = 0.0  # Für gewichtetes Mittel
            
            for i in range(len(obs_bv_dereddened)):
                # Finde nächsten Punkt auf der Isochrone in B-V
                color_diff = np.abs(iso_bv_trimmed - obs_bv_dereddened[i])
                if len(color_diff) == 0:
                    continue
                    
                idx = np.argmin(color_diff)
                
                # Nur verwenden wenn Farbe nahe genug (innerhalb 0.15 mag in B-V)
                if color_diff[idx] < 0.15:
                    v_diff = obs_v_dereddened[i] - iso_v_shifted[idx]
                    
                    # NUR Sterne verwenden, die am BLAUEN RAND liegen
                    # v_diff < 0: Stern ist HELLER als Isochrone (physikalisch unrealistisch für ZAMS)
                    # v_diff > 0: Stern ist SCHWÄCHER als Isochrone (Binäre, Post-MS, etc.)
                    # Wir wollen nur Sterne nahe der ZAMS (kleine positive v_diff)
                    if -0.5 <= v_diff <= 0.5:  # Stern ist max 0.5 mag vom Isochrone entfernt
                        # Gewichte Sterne nahe der Isochrone stärker
                        # Je näher am blauen Rand (kleine v_diff), desto höher das Gewicht
                        weight = np.exp(-abs(v_diff) / 0.3)
                        chi2 += weight * v_diff**2
                        total_weight += weight
                        n_matched += 1
            
            if n_matched > 5 and total_weight > 0:  # Mindestens 5 gematchte Sterne
                chi2 /= total_weight  # Gewichtetes mittleres Chi² (robuster!)
                
                # WICHTIG: Regularisierung gegen unphysikalische Entfernungen
                # Bevorzuge Entfernungen im Bereich 800-2000 pc (typisch für offene Haufen)
                # Prior: Gaussian mit Zentrum bei 1300 pc, Breite 500 pc (engerer Prior!)
                distance_pc = 10**((dm+5)/5)
                distance_prior = np.exp(-((distance_pc - 1300)**2) / (2 * 500**2))
                
                # Kombiniere Chi² mit Prior (sehr starker Prior, Faktor 1.5)
                # Dies bevorzugt physikalisch plausible Entfernungen sehr stark
                # und verhindert Sprünge zu unphysikalischen lokalen Minima
                chi2_with_prior = chi2 * (1 + 1.5 * (1 - distance_prior))
                
                # Speichere für Debug
                if chi2_curve is not None:
                    chi2_curve.append(chi2_with_prior)
                    dm_curve.append(dm)
                
                if chi2_with_prior < best_chi2:
                    best_chi2 = chi2_with_prior
                    best_dm = dm
        
        # Debug: Zeige Chi² Kurve bei fixiertem Alter
        if chi2_curve is not None and len(chi2_curve) > 0:
            chi2_arr = np.array(chi2_curve)
            dm_arr = np.array(dm_curve)
            
            # Finde alle lokalen Minima
            local_minima = []
            for i in range(1, len(chi2_arr) - 1):
                if chi2_arr[i] < chi2_arr[i-1] and chi2_arr[i] < chi2_arr[i+1]:
                    local_minima.append((dm_arr[i], chi2_arr[i], 10**((dm_arr[i]+5)/5)))
            
            print(f"\n  === Chi² Analyse für log(age)={age_log:.2f} ===")
            print(f"  Gefundene lokale Minima: {len(local_minima)}")
            for dm_min, chi2_min, dist_min in local_minima:
                print(f"    DM={dm_min:.2f}, d={dist_min:.0f} pc, χ²/N={chi2_min:.4f}")
            print(f"  Globales Minimum: DM={best_dm:.2f}, d={10**((best_dm+5)/5):.0f} pc, χ²/N={best_chi2:.4f}\n")
        
        # Debug: Zeige beste Anpassung für dieses Alter
        if best_dm is not np.nan and not np.isnan(best_dm):
            print(f"  log(age)={age_log:.2f}: DM={best_dm:.2f}, d={10**((best_dm+5)/5):.0f} pc, χ²/N={best_chi2:.3f}")
        
        chi2_values.append(best_chi2)
        distance_moduli.append(best_dm)
    
    # Finde bestes Alter
    chi2_values = np.array(chi2_values)
    distance_moduli = np.array(distance_moduli)
    
    valid = np.isfinite(chi2_values)
    if not np.any(valid):
        print("FEHLER: Keine gültige Anpassung gefunden!")
        return None
    
    best_idx = np.argmin(chi2_values[valid])
    best_age_log = ages_log[valid][best_idx]
    best_dm = distance_moduli[valid][best_idx]
    best_chi2 = chi2_values[valid][best_idx]
    
    # Lade beste Isochrone nochmal
    iso_bv_best, iso_v_best = download_isochrone_cmd(best_age_log)
    
    # Berechne Entfernung in Parsec
    # m - M = 5 * log10(d) - 5
    # d = 10^((m-M+5)/5)
    distance_pc = 10**((best_dm + 5) / 5.0)
    distance_ly = distance_pc * 3.26156  # Lichtjahre
    
    print(f"\n=== Ergebnisse ===")
    if fixed_age_myr is not None:
        print(f"Alter (fixiert): {fixed_age_myr} Millionen Jahre (log(age) = {best_age_log:.3f})")
    else:
        print(f"Bestes Alter: log(age) = {best_age_log:.2f} ({10**best_age_log / 1e6:.1f} Millionen Jahre)")
    print(f"Entfernungsmodul: (m-M)_0 = {best_dm:.2f} mag")
    print(f"Entfernung: d = {distance_pc:.0f} pc = {distance_pc/1000:.2f} kpc = {distance_ly:.0f} Lichtjahre")
    print(f"Chi²/N = {best_chi2:.3f}")
    print(f"\nZum Vergleich:")
    print(f"  1000 pc → DM = {5*np.log10(1000)-5:.2f} mag")
    print(f"  1300 pc → DM = {5*np.log10(1300)-5:.2f} mag")
    print(f"  2000 pc → DM = {5*np.log10(2000)-5:.2f} mag")
    
    # Plotte Chi² vs. Alter - IMMER erstellen und speichern
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7))
    
    # Chi² Plot
    ax1.plot(ages_log[valid], chi2_values[valid], 'b.-')
    ax1.axvline(best_age_log, color='r', linestyle='--', label=f'Best: log(age)={best_age_log:.2f}')
    ax1.set_xlabel('log(Alter / Jahre)')
    ax1.set_ylabel('χ²/N')
    ax1.set_title('Isochronen-Anpassung: Chi² vs. Alter')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # CMD mit bester Isochrone - zeige VOLLE Isochrone im Plot
    iso_bv_best, iso_v_best = download_isochrone_cmd(best_age_log, main_sequence_only=False)
    
    ax2.scatter(obs_bv_dereddened, obs_v_dereddened, s=20, alpha=0.6, label='Beobachtungen (dereddened)')
    if iso_bv_best is not None:
        # Bestimme Bereich wo Sterne sind für Achsenbegrenzung
        bv_min, bv_max = np.min(obs_bv_dereddened), np.max(obs_bv_dereddened)
        v_min, v_max = np.min(obs_v_dereddened), np.max(obs_v_dereddened)
        
        # Erweitere um 10% für bessere Darstellung
        bv_margin = (bv_max - bv_min) * 0.1
        v_margin = (v_max - v_min) * 0.1
        
        # Filter Isochrone auf sichtbaren Bereich (etwas erweitert)
        iso_mask = ((iso_bv_best >= bv_min - bv_margin) & 
                    (iso_bv_best <= bv_max + bv_margin) &
                    (iso_v_best + best_dm >= v_min - v_margin) &
                    (iso_v_best + best_dm <= v_max + v_margin))
        
        ax2.plot(iso_bv_best[iso_mask], (iso_v_best + best_dm)[iso_mask], 
                'r-', linewidth=2, label=f'Isochrone (log(age)={best_age_log:.2f}, DM={best_dm:.2f})')
        
        # Setze Achsenbegrenzung auf Sternbereich
        ax2.set_xlim(bv_min - bv_margin, bv_max + bv_margin)
        ax2.set_ylim(v_max + v_margin, v_min - v_margin)  # Invertiert!
    ax2.set_xlabel('(B-V)₀ [mag]')
    ax2.set_ylabel('V₀ [mag]')
    ax2.set_title(f'CMD mit bester Isochrone\nd = {distance_pc:.0f} pc, Alter = {10**best_age_log/1e6:.1f} Myr')
    ax2.legend()
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    
    # Speichere immer
    if output_folder:
        fit_plot_path = get_output_path(output_folder, 'isochrone_fit.png')
        fig.savefig(fit_plot_path, dpi=150)
        print(f"Anpassungs-Plot gespeichert: {fit_plot_path}")
    
    plt.show()
    
    # Gebe Dictionary mit allen Parametern zurück
    return {
        'age_log': best_age_log,
        'dm': best_dm,
        'distance_pc': distance_pc,
        'chi2': best_chi2,
        'iso_bv': iso_bv_best,
        'iso_v': iso_v_best,
        'ebv': ebv,
        'obs_bv_dereddened': obs_bv_dereddened,
        'obs_v_dereddened': obs_v_dereddened
    }


def plot_color_magnitude(matched, tabB, tabV, paths, n_total_B: int = None, n_total_V: int = None, 
                         calibrated: bool = False, mags_B=None, mags_V=None, colors=None,
                         save: Optional[str] = None, 
                         isochrone_params: Optional[dict] = None,
                         outlier_sigma: float = 3.0):
    # Verwende entweder kalibrierte oder instrumentale Magnituden
    if calibrated and mags_B is not None and mags_V is not None and colors is not None:
        magsB = mags_B
        magsV = mags_V
        color = colors
        ylabel = 'V (mag)'
        title_prefix = 'Kalibriertes Color-Magnitude'
    else:
        # Verwende hintergrundkorrigierte Flüsse für instrumentale Photometrie
        sumsB = np.asarray(tabB['aperture_sum_bkgsub'])
        sumsV = np.asarray(tabV['aperture_sum_bkgsub'])
        magsB = []
        magsV = []
        for iB, iV in matched:
            magB = instrumental_magnitude(sumsB[iB])
            magV = instrumental_magnitude(sumsV[iV])
            if magB is None or magV is None:
                continue
            magsB.append(magB)
            magsV.append(magV)
        magsB = np.array(magsB)
        magsV = np.array(magsV)
        if len(magsB) == 0:
            print("Keine validen Matches mit positivem Fluss gefunden.")
            return
        color = magsB - magsV  # B - V
        ylabel = 'V (instrumental mag)'
        title_prefix = 'Color-Magnitude'
    
    # Entferne Ausreißer basierend auf Sigma-Clipping für Farbe
    # Realistische B-V Werte liegen zwischen etwa -0.3 und +2.0
    color_median = np.median(color)
    color_std = np.std(color)
    color_mask = np.abs(color - color_median) <= outlier_sigma * color_std
    
    n_outliers = np.sum(~color_mask)
    if n_outliers > 0:
        print(f"\n⚠️  Entferne {n_outliers} Farb-Ausreißer (>{outlier_sigma:.1f}σ)")
        print(f"    Original: B-V = {np.min(color):.2f}...{np.max(color):.2f}")
        magsB = magsB[color_mask]
        magsV = magsV[color_mask]
        color = color[color_mask]
        print(f"    Nach Filter: B-V = {np.min(color):.2f}...{np.max(color):.2f}")
    
    # Debug-Ausgabe für Diagnose
    print(f"\n=== FHD Statistik ===")
    print(f"Anzahl Sterne: {len(color)}")
    print(f"V-Magnitude: min={np.min(magsV):.2f}, max={np.max(magsV):.2f}, median={np.median(magsV):.2f}")
    print(f"B-V Farbe: min={np.min(color):.2f}, max={np.max(color):.2f}, median={np.median(color):.2f}")
    
    # Wenn unkalibriert und B-V sehr negativ: warne über Nullpunkt-Problem
    if not calibrated and np.median(color) < -0.5:
        offset = np.median(color)
        print(f"\n⚠️  WARNUNG: B-V Median = {offset:.2f} ist sehr negativ!")
        print(f"    Dies deutet auf unterschiedliche Nullpunkte zwischen B und V hin.")
        print(f"    Für kalibrierte Magnituden verwenden Sie: --calibrate")
        print(f"    Oder testen Sie grob: Die Hauptreihe hätte B-V ≈ {offset + 1.9:.2f} wenn")
        print(f"    Sie {abs(offset):.2f} mag vom B-V Wert abziehen würden.")
    
    print(f"\nErwartung: Helle Sterne (kleine V) oben, schwache Sterne (große V) unten")
    print(f"Erwartung: Blaue Sterne (negative/kleine B-V) links, rote Sterne (positive/große B-V) rechts")

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(color, magsV, s=30, c='tab:blue', alpha=0.8, label='Beobachtungen')
    
    # Füge Isochrone hinzu wenn Parameter vorhanden
    if isochrone_params is not None and calibrated:
        iso_bv = isochrone_params['iso_bv']
        iso_v = isochrone_params['iso_v']
        dm = isochrone_params['dm']
        age_log = isochrone_params['age_log']
        ebv = isochrone_params['ebv']
        distance_pc = isochrone_params['distance_pc']
        
        # WICHTIG: Wenn Isochronen-Parameter vorhanden, sind die Daten bereits dereddened!
        # Isochrone ist intrinsisch, also direkt vergleichbar
        iso_v_shifted = iso_v + dm  # Nur Entfernungsmodul, KEINE Extinktion
        
        # Filter Isochrone auf sichtbaren Bereich
        bv_min, bv_max = np.min(color), np.max(color)
        v_min, v_max = np.min(magsV), np.max(magsV)
        bv_margin = (bv_max - bv_min) * 0.15
        v_margin = (v_max - v_min) * 0.15
        
        iso_mask = ((iso_bv >= bv_min - bv_margin) & 
                    (iso_bv <= bv_max + bv_margin) &
                    (iso_v_shifted >= v_min - v_margin) &
                    (iso_v_shifted <= v_max + v_margin))
        
        if np.any(iso_mask):
            ax.plot(iso_bv[iso_mask], iso_v_shifted[iso_mask], 
                   'r-', linewidth=2.5, alpha=0.8,
                   label=f'Isochrone ({10**age_log/1e6:.1f} Myr)')
        
        # Aktualisiere Titel mit physikalischen Parametern
        title = f'{title_prefix} (dereddened)\n'
        title += f'd = {distance_pc:.0f} pc, Alter = {10**age_log/1e6:.1f} Myr, E(B-V) = {ebv:.3f}'
    else:
        # Titel mit Info über gefilterte Sterne
        title = f'{title_prefix} ({os.path.basename(paths[0])} vs {os.path.basename(paths[1])})'
        if n_total_B is not None and n_total_V is not None:
            title += f'\n{len(color)} Sterne (gefiltert aus {n_total_B}/{n_total_V} Gesamtquellen)'
    
    ax.set_xlabel('(B-V)₀ [mag]')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    # Invertiere y-Achse für astronomische Orientierung:
    # - y-Achse (V): kleine Magnituden (helle Sterne) oben
    # - x-Achse (B-V): negative Werte (blaue Sterne) links (Standard)
    ax.invert_yaxis()  # Heller oben
    
    # Füge ein Raster hinzu und markiere Null-Linie für B-V
    ax.axvline(x=0, color='gray', linestyle=':', alpha=0.5)
    ax.grid(alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    if save:
        try:
            fig.savefig(save, dpi=150)
            print(f"Color-Magnitude Abbildung gespeichert: {save}")
        except Exception as e:
            print(f"Fehler beim Speichern des Color-Magnitude Diagramms: {e}", file=sys.stderr)
    plt.show()


def compute_instrumental_cmd_points(matched, tabB, tabV, outlier_sigma: float = 3.0):
    """Berechne instrumentale CMD-Punkte (B-V, V) inkl. Sigma-Clipping der Farben."""
    if matched is None or tabB is None or tabV is None or len(matched) == 0:
        return None

    sumsB = np.asarray(tabB['aperture_sum_bkgsub'])
    sumsV = np.asarray(tabV['aperture_sum_bkgsub'])
    magsB = []
    magsV = []
    x_coords = []
    y_coords = []

    for iB, iV in matched:
        magB = instrumental_magnitude(sumsB[iB])
        magV = instrumental_magnitude(sumsV[iV])
        if magB is None or magV is None:
            continue
        magsB.append(magB)
        magsV.append(magV)
        x_coords.append(tabB['x'][iB])
        y_coords.append(tabB['y'][iB])

    magsB = np.array(magsB)
    magsV = np.array(magsV)
    x_coords = np.array(x_coords)
    y_coords = np.array(y_coords)
    if len(magsB) == 0:
        return None

    color = magsB - magsV
    if len(color) > 1:
        color_median = np.median(color)
        color_std = np.std(color)
        if color_std > 0:
            color_mask = np.abs(color - color_median) <= outlier_sigma * color_std
        else:
            color_mask = np.ones_like(color, dtype=bool)
    else:
        color_mask = np.ones_like(color, dtype=bool)

    return {
        'color': color[color_mask],
        'mag_v': magsV[color_mask],
        'x': x_coords[color_mask],
        'y': y_coords[color_mask],
    }


def find_nearest_source_index(tab, x_ref: float, y_ref: float, max_distance: Optional[float] = None):
    """Finde den naechsten Stern in tab zu einer Referenzposition (x_ref, y_ref)."""
    if tab is None or len(tab) == 0:
        return None

    coords = np.vstack((np.array(tab['x']), np.array(tab['y']))).T
    ref = np.array([x_ref, y_ref])
    dists = np.linalg.norm(coords - ref, axis=1)
    idx = int(np.argmin(dists))

    if max_distance is not None and dists[idx] > max_distance:
        return None
    return idx


def compute_calibrated_cmd_points(matched, tabB, tabV, ref_index_B: int,
                                  mag_B_true: float, mag_V_true: float,
                                  outlier_sigma: float = 3.0,
                                  ebv: Optional[float] = None):
    """Berechne kalibrierte CMD-Punkte; optional mit Dereddening via E(B-V)."""
    calib = calibrate_magnitudes(matched, tabB, tabV, ref_index_B, mag_B_true, mag_V_true)
    if calib is None:
        return None

    magsB, magsV, colors, valid_matches = calib
    if magsB is None or len(magsB) == 0:
        return None

    x_coords = np.array([tabB['x'][iB] for iB, _ in valid_matches])
    y_coords = np.array([tabB['y'][iB] for iB, _ in valid_matches])

    if len(colors) > 1:
        color_median = np.median(colors)
        color_std = np.std(colors)
        if color_std > 0:
            color_mask = np.abs(colors - color_median) <= outlier_sigma * color_std
        else:
            color_mask = np.ones_like(colors, dtype=bool)
    else:
        color_mask = np.ones_like(colors, dtype=bool)

    color_plot = colors[color_mask]
    magsV_plot = magsV[color_mask]

    if ebv is not None:
        color_plot = color_plot - ebv
        magsV_plot = magsV_plot - 3.1 * ebv

    return {
        'color': color_plot,
        'mag_v': magsV_plot,
        'x': x_coords[color_mask],
        'y': y_coords[color_mask],
    }


def plot_cmd_method_comparison(aperture_points, psf_points, paths, save: Optional[str] = None,
                               calibrated: bool = False, dereddened: bool = False):
    """Plotte CMD-Punkte fuer aperture und psf in verschiedenen Farben."""
    if aperture_points is None and psf_points is None:
        print("Keine CMD-Daten fuer den Methodenvergleich verfuegbar.")
        return

    fig, ax = plt.subplots(figsize=(7, 6))
    plotted_any = False
    n_connections = 0
    delta_color_values = []
    delta_mag_values = []

    if aperture_points is not None:
        color_ap = aperture_points['color']
        magsV_ap = aperture_points['mag_v']
        if len(color_ap) > 0:
            ax.scatter(color_ap, magsV_ap, s=28, c='tab:blue', alpha=0.7,
                       label=f'aperture ({len(color_ap)} Sterne)')
            plotted_any = True

    if psf_points is not None:
        color_psf = psf_points['color']
        magsV_psf = psf_points['mag_v']
        if len(color_psf) > 0:
            ax.scatter(color_psf, magsV_psf, s=28, c='tab:orange', alpha=0.7,
                       label=f'psf ({len(color_psf)} Sterne)')
            plotted_any = True

    if aperture_points is not None and psf_points is not None:
        ap_xy = np.column_stack((aperture_points['x'], aperture_points['y']))
        psf_xy = np.column_stack((psf_points['x'], psf_points['y']))
        used_psf = set()
        max_position_delta = 2.0

        for i_ap, xy_ap in enumerate(ap_xy):
            if len(psf_xy) == 0:
                break
            dists = np.linalg.norm(psf_xy - xy_ap, axis=1)
            i_psf = int(np.argmin(dists))
            if dists[i_psf] > max_position_delta or i_psf in used_psf:
                continue

            used_psf.add(i_psf)
            delta_color_values.append(psf_points['color'][i_psf] - aperture_points['color'][i_ap])
            delta_mag_values.append(psf_points['mag_v'][i_psf] - aperture_points['mag_v'][i_ap])
            ax.plot(
                [aperture_points['color'][i_ap], psf_points['color'][i_psf]],
                [aperture_points['mag_v'][i_ap], psf_points['mag_v'][i_psf]],
                color='0.55', alpha=0.45, linewidth=1.6, zorder=1
            )
            n_connections += 1

    if not plotted_any:
        print("Keine validen CMD-Punkte fuer den Methodenvergleich gefunden.")
        plt.close(fig)
        return

    if n_connections > 0:
        print(f"  Verbindungssegmente fuer {n_connections} gemeinsame Sterne gezeichnet.")
        delta_color_values = np.array(delta_color_values)
        delta_mag_values = np.array(delta_mag_values)

        rms_color = float(np.sqrt(np.mean(delta_color_values**2)))
        rms_mag = float(np.sqrt(np.mean(delta_mag_values**2)))
        mean_color = float(np.mean(delta_color_values))
        mean_mag = float(np.mean(delta_mag_values))

        color_tendency = 'PSF tendenziell roeter' if mean_color > 0 else 'PSF tendenziell blauer'
        if np.isclose(mean_color, 0.0, atol=1e-4):
            color_tendency = 'keine klare Farb-Tendenz'

        mag_tendency = 'PSF tendenziell dunkler' if mean_mag > 0 else 'PSF tendenziell heller'
        if np.isclose(mean_mag, 0.0, atol=1e-4):
            mag_tendency = 'keine klare Helligkeits-Tendenz'

        print(f"  RMS Delta(B-V) [psf-aperture] = {rms_color:.4f} mag")
        print(f"  RMS Delta(V) [psf-aperture] = {rms_mag:.4f} mag")
        print(f"  Mittel Delta(B-V) = {mean_color:+.4f} mag -> {color_tendency}")
        print(f"  Mittel Delta(V) = {mean_mag:+.4f} mag -> {mag_tendency}")

        stats_text = (
            f'Gemeinsame Sterne: {n_connections}\n'
            f'RMS Delta(B-V): {rms_color:.4f} mag\n'
            f'RMS Delta(V): {rms_mag:.4f} mag\n'
            f'Mittel Delta(B-V): {mean_color:+.4f} mag\n'
            f'{color_tendency}\n'
            f'Mittel Delta(V): {mean_mag:+.4f} mag\n'
            f'{mag_tendency}'
        )
        ax.text(
            0.02, 0.98, stats_text,
            transform=ax.transAxes,
            va='top', ha='left', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, edgecolor='0.7')
        )

    if calibrated:
        if dereddened:
            xlabel = '(B-V)0 [mag]'
            ylabel = 'V0 [mag]'
            title_mode = 'kalibriert + dereddened'
        else:
            xlabel = '(B-V) [mag]'
            ylabel = 'V [mag]'
            title_mode = 'kalibriert'
    else:
        xlabel = '(B-V) [mag]'
        ylabel = 'V (instrumental mag)'
        title_mode = 'instrumentell'

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(
        f'CMD-Vergleich aperture vs psf ({title_mode})\n'
        f'{os.path.basename(paths[0])} vs {os.path.basename(paths[1])}'
    )
    ax.axvline(x=0, color='gray', linestyle=':', alpha=0.5)
    ax.grid(alpha=0.3)
    ax.legend()
    ax.invert_yaxis()

    plt.tight_layout()
    if save:
        try:
            fig.savefig(save, dpi=150)
            print(f"CMD-Vergleich gespeichert: {save}")
        except Exception as e:
            print(f"Fehler beim Speichern des CMD-Vergleichs: {e}", file=sys.stderr)
    plt.show()


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description='Zeige zwei FITS-Dateien aus einem Ordner an und führe Aperturphotometrie durch')
    parser.add_argument('folder', help='Pfad zum Ordner mit FITS-Dateien')
    parser.add_argument('--indices', '-i', type=int, nargs='*', help='Index(e) der anzuzeigenden Dateien (0-basiert). Wenn nicht gesetzt, werden die ersten zwei verwendet.')
    parser.add_argument('--save', '-s', help='Optional: Pfad (z.B. out.png) um die erzeugte Abbildung der Frames zu speichern')
    parser.add_argument('--save-cm', help='Optional: Pfad (z.B. cm.png) um das Color-Magnitude-Diagramm zu speichern')
    parser.add_argument('--compare-cmd-methods', action='store_true',
                        help='Vergleiche CMDs von aperture und psf in einem gemeinsamen Plot mit unterschiedlichen Farben.')
    parser.add_argument('--save-cmd-compare', type=str, default=None,
                        help='Optional: Pfad (z.B. cmd_compare.png) zum Speichern des CMD-Methodenvergleichs.')
    parser.add_argument('--verbose', '-v', action='store_true', help='Mehr Ausgaben')
    parser.add_argument('--aperture', type=float, default=None, 
                       help='Aperturradius in Pixel. Wenn nicht angegeben: automatisch 2.5×FWHM (empfohlen!)')
    parser.add_argument('--fwhm', type=float, default=3.0, help='FWHM in Pixel für Quellenfindung (Standard 3.0)')
    parser.add_argument('--threshold', type=float, default=5.0, help='Schwelle in sigma für Quellenfindung (Standard 5.0)')
    parser.add_argument('--photometry-method', type=str, choices=['aperture', 'psf'], default='aperture',
                        help='Photometrie-Methode: aperture (Standard) oder psf (besser für dichte Sternfelder)')
    parser.add_argument('--maxsep', type=float, default=5.0, help='Maximale Trennung in Pixel zum Matchen (Standard 5.0)')
    parser.add_argument('--filters', nargs=2, default=['B', 'V'], help='Labels der beiden Bilder (z.B. B V). Erstes Bild -> erstes Label.')
    parser.add_argument('--force-filters', action='store_true',
                        help='Erzwinge --filters auch wenn FITS-Header Filter-Info haben (überschreibt Header)')
    parser.add_argument('--cluster-radius', type=float, default=None, 
                        help='Radius in Pixel um das Zentrum für Sternhaufen-Filterung. Standard: 40%% der y-Bildgröße. 0 = deaktiviert.')
    parser.add_argument('--cluster-center', type=float, nargs=2, default=None, metavar=('X', 'Y'),
                        help='Zentrum (x, y) in Pixel für Sternhaufen-Filterung. Standard: Bildmitte.')
    parser.add_argument('--annulus-radii', type=float, nargs=2, default=None, metavar=('R_INNER', 'R_OUTER'),
                        help=('Nur Sterne in einem Ring (Annulus) zwischen innerem (R_INNER) und äußerem '
                              '(R_OUTER) Radius in Pixel um das Zentrum berücksichtigen. '
                              'Überschreibt --cluster-radius.'))
    parser.add_argument('--calibrate', action='store_true',
                        help='Aktiviere interaktive Kalibrierung mit Referenzstern')
    parser.add_argument('--calibrate-file', type=str, default=None,
                        help='Pfad zur Kalibrierungsdatei (zum Speichern bei --calibrate oder Laden)')
    parser.add_argument('--fit-isochrone', action='store_true',
                        help='Bestimme Entfernung und Alter durch Isochronen-Anpassung (benötigt --calibrate)')
    parser.add_argument('--ebv', type=float, default=None,
                        help='Rötung E(B-V) in mag für Isochronen-Anpassung (wird abgefragt falls nicht angegeben)')
    parser.add_argument('--age-range', type=float, nargs=2, default=[6.8, 7.6],
                        help='Altersbereich für Isochronen-Anpassung als log(age/Jahre) (Standard: 6.8 7.6 = 6-40 Myr)')
    parser.add_argument('--age', type=float, default=None,
                        help='Fixiertes Alter in Millionen Jahren (z.B. 10.6 für NGC869, 25 für M36). Wenn gesetzt, wird nur die Entfernung bestimmt.')
    parser.add_argument('--outlier-sigma', type=float, default=3.0,
                        help='Sigma-Schwelle für Ausreißer-Filterung in B-V Farbe (Standard: 3.0). Höhere Werte = weniger Filterung.')
    args = parser.parse_args(argv)

    folder = args.folder
    if not os.path.isdir(folder):
        print(f"Ordner existiert nicht: {folder}", file=sys.stderr)
        return 2
    
    # Bestimme Ausgabepfade: Falls relative Pfade oder nur Dateinamen angegeben,
    # speichere im FITS-Ordner. Absolute Pfade bleiben unverändert.
    save_images_path = None
    if args.save:
        if os.path.isabs(args.save):
            save_images_path = args.save
        else:
            save_images_path = get_output_path(folder, args.save)
    
    save_cmd_path = None
    if args.save_cm:
        if os.path.isabs(args.save_cm):
            save_cmd_path = args.save_cm
        else:
            save_cmd_path = get_output_path(folder, args.save_cm)

    save_cmd_compare_path = None
    if args.save_cmd_compare:
        if os.path.isabs(args.save_cmd_compare):
            save_cmd_compare_path = args.save_cmd_compare
        else:
            save_cmd_compare_path = get_output_path(folder, args.save_cmd_compare)
    
    calibrate_file_path = None
    if args.calibrate_file:
        if os.path.isabs(args.calibrate_file):
            calibrate_file_path = args.calibrate_file
        else:
            calibrate_file_path = get_output_path(folder, args.calibrate_file)
    elif args.calibrate:
        # Kein expliziter Pfad angegeben, verwende Standard im FITS-Ordner
        calibrate_file_path = get_output_path(folder, 'calib.txt')
    
    if args.verbose:
        if save_images_path:
            print(f"Bilder werden gespeichert in: {save_images_path}")
        if save_cmd_path:
            print(f"CMD wird gespeichert in: {save_cmd_path}")
        if save_cmd_compare_path:
            print(f"CMD-Methodenvergleich wird gespeichert in: {save_cmd_compare_path}")
        if calibrate_file_path:
            print(f"Kalibrierungsdatei: {calibrate_file_path}")

    files = find_fits_files(folder)
    if args.verbose:
        print(f"Gefundene FITS-Dateien ({len(files)}):")
        for i, f in enumerate(files):
            print(f"  {i}: {f}")

    if len(files) < 2 and not args.indices:
        print(f"Es wurden nur {len(files)} FITS-Datei(en) gefunden. Mindestens 2 werden erwartet.", file=sys.stderr)
        if len(files) == 0:
            return 1

    if args.indices:
        selected = []
        for idx in args.indices:
            if idx < 0 or idx >= len(files):
                print(f"Index {idx} außerhalb der gefundenen Dateien.", file=sys.stderr)
                return 1
            selected.append(files[idx])
    else:
        if len(files) == 0:
            print("Keine FITS-Dateien gefunden.", file=sys.stderr)
            return 1
        # Standard: die ersten zwei
        selected = files[:2]
        if len(files) > 2 and args.verbose:
            print(f"Mehr als 2 FITS-Dateien gefunden; zeige die ersten beiden: {os.path.basename(selected[0])}, {os.path.basename(selected[1])}")

    # Lade Bilder
    images = []
    datas = []
    headers = []
    for p in selected:
        try:
            data, header = load_fits_image(p)
        except Exception as e:
            print(str(e), file=sys.stderr)
            data, header = None, None
        images.append((data, header))
        if data is not None:
            datas.append(data)
        else:
            print(f"Überspringe Photometrie für {p}: keine gültigen Daten")
            datas.append(None)
        headers.append(header)

    # Bestimme Zentrum und Radius für Sternhaufen-Filterung (vor der Verwendung!)
    cluster_center = None
    cluster_radius = None
    
    # Berechne cluster_radius und cluster_center basierend auf den Argumenten
    if args.cluster_radius is not None and args.cluster_radius > 0:
        cluster_radius = args.cluster_radius
    elif args.cluster_radius != 0:  # None oder nicht explizit 0
        # Standard: 40% der y-Bildgröße
        if datas[0] is not None:
            cluster_radius = datas[0].shape[0] * 0.4
        elif datas[1] is not None:
            cluster_radius = datas[1].shape[0] * 0.4
    
    if cluster_radius is not None and cluster_radius > 0:
        if args.cluster_center is not None:
            cluster_center = tuple(args.cluster_center)
        else:
            # Standard: Bildmitte
            if datas[0] is not None:
                cluster_center = (datas[0].shape[1] / 2.0, datas[0].shape[0] / 2.0)
            elif datas[1] is not None:
                cluster_center = (datas[1].shape[1] / 2.0, datas[1].shape[0] / 2.0)
        
        if args.verbose and cluster_center is not None:
            print(f"Sternhaufen-Filterung: Radius={cluster_radius:.1f} Pixel, Zentrum=({cluster_center[0]:.1f}, {cluster_center[1]:.1f})")
    
    # Annulus-Filterung verarbeiten (überschreibt ggf. cluster_radius)
    annulus_inner = None
    if args.annulus_radii is not None:
        r_inner, r_outer = args.annulus_radii
        if r_inner < 0 or r_outer <= 0:
            print("FEHLER: Annulus-Radien müssen positiv sein.", file=sys.stderr)
            return 1
        if r_inner >= r_outer:
            print(f"FEHLER: Innerer Radius ({r_inner:.1f}) muss kleiner als äußerer Radius ({r_outer:.1f}) sein.", file=sys.stderr)
            return 1
        annulus_inner = r_inner
        cluster_radius = r_outer
        # Stelle sicher, dass cluster_center gesetzt ist
        if cluster_center is None:
            if datas[0] is not None:
                cluster_center = (datas[0].shape[1] / 2.0, datas[0].shape[0] / 2.0)
            elif datas[1] is not None:
                cluster_center = (datas[1].shape[1] / 2.0, datas[1].shape[0] / 2.0)
        print(f"Annulus-Filterung: R_inner={annulus_inner:.1f}, R_outer={cluster_radius:.1f} Pixel, "
              f"Zentrum=({cluster_center[0]:.1f}, {cluster_center[1]:.1f})")
    
    # Photometrie durchführen
    try:
        method = args.photometry_method if hasattr(args, 'photometry_method') else 'aperture'
        if args.verbose:
            print(f"Verwende {method}-Photometrie...")
        
        result1 = detect_and_photometry_wrapper(datas[0], method=method, fwhm=args.fwhm, 
                                               threshold_sigma=args.threshold, 
                                               aperture_radius=args.aperture) if datas[0] is not None else (None, None)
        result2 = detect_and_photometry_wrapper(datas[1], method=method, fwhm=args.fwhm, 
                                               threshold_sigma=args.threshold, 
                                               aperture_radius=args.aperture) if datas[1] is not None else (None, None)
        
        # Entpacke Ergebnisse
        tab1, aperture_used_1 = result1
        tab2, aperture_used_2 = result2
        
        # Verwende den ersten nicht-None Aperturwert für die Plots
        aperture_for_plot = aperture_used_1 if aperture_used_1 is not None else aperture_used_2
        
        if args.verbose and (tab1 is not None or tab2 is not None):
            n1 = len(tab1) if tab1 is not None else 0
            n2 = len(tab2) if tab2 is not None else 0
            print(f"Gefundene Quellen: {n1} im ersten, {n2} im zweiten Bild")
    except Exception as e:
        print(f"Fehler bei der Photometrie: {e}", file=sys.stderr)
        return 1
    
    if tab1 is None or tab2 is None:
        print("In mindestens einem Bild wurden keine Quellen gefunden.")
        return 1

    # Führe Matching durch um den Offset zu bekommen
    matched, image_offset = match_sources(tab1, tab2, max_sep=args.maxsep, 
                           filter_radius=cluster_radius, center=cluster_center,
                           annulus_inner=annulus_inner)
    if args.verbose:
        print(f"Gefundene Matches: {len(matched)}")

    # Zeige die Frames mit Aperturen (NACH Matching, damit Offset bekannt ist)
    # Übergebe matched_indices damit nur gematchte Sterne angezeigt werden
    plot_images(images, selected, phot_tables=[tab1, tab2], aperture_radius=aperture_for_plot, 
                cluster_radius=cluster_radius, cluster_center=cluster_center, 
                image_offset=image_offset, matched_indices=matched, annulus_inner=annulus_inner,
                save=save_images_path)

    # Bestimme Filter aus FITS-Headern
    filter1_from_header = get_filter_from_header(headers[0]) if len(headers) > 0 else None
    filter2_from_header = get_filter_from_header(headers[1]) if len(headers) > 1 else None
    
    print(f"\n=== Filter-Zuordnung ===")
    print(f"Erstes Bild: {os.path.basename(selected[0])}")
    if filter1_from_header:
        print(f"  FITS-Header: FILTER = '{filter1_from_header}'")
    else:
        print(f"  FITS-Header: Kein Filter-Keyword gefunden")
    
    print(f"Zweites Bild: {os.path.basename(selected[1])}")
    if filter2_from_header:
        print(f"  FITS-Header: FILTER = '{filter2_from_header}'")
    else:
        print(f"  FITS-Header: Kein Filter-Keyword gefunden")
    
    # Bestimme Filter-Zuordnung
    # Priorität 1: --force-filters überschreibt alles
    # Priorität 2: FITS-Header
    # Priorität 3: Command-line --filters Argumente
    
    filters_swapped = False  # Track ob wir vertauscht haben
    
    if args.force_filters:
        # Erzwinge --filters auch wenn Header vorhanden
        label1, label2 = args.filters
        print(f"\n⚠ --force-filters aktiv: Ignoriere FITS-Header")
        print(f"Verwende --filters Argumente: '{label1}' und '{label2}'")
        
        if label1.upper() == 'B' and label2.upper() == 'V':
            tabB = tab1
            tabV = tab2
            print(f"  → Zuordnung: Bild 1 = B, Bild 2 = V")
        elif label1.upper() == 'V' and label2.upper() == 'B':
            tabB = tab2
            tabV = tab1
            filters_swapped = True
            print(f"  → Zuordnung: Bild 1 = V, Bild 2 = B (vertauscht)")
        else:
            print(f"  Unklare Labels. Verwende Flux-basierte Heuristik...")
            tabB, tabV = assign_filters_by_flux(tab1, tab2)
    elif filter1_from_header and filter2_from_header:
        # Verwende Header-Informationen
        if filter1_from_header == 'B' and filter2_from_header == 'V':
            tabB = tab1
            tabV = tab2
            print(f"\n✓ Filter aus FITS-Header: Bild 1 = B, Bild 2 = V")
        elif filter1_from_header == 'V' and filter2_from_header == 'B':
            tabB = tab2
            tabV = tab1
            filters_swapped = True
            print(f"\n✓ Filter aus FITS-Header: Bild 1 = V, Bild 2 = B")
            print(f"  → Vertausche Zuordnung: tabB = Bild 2, tabV = Bild 1")
        else:
            print(f"\n⚠ Unerwartete Filter in Header: '{filter1_from_header}' und '{filter2_from_header}'")
            print(f"  Verwende Flux-basierte Heuristik...")
            tabB, tabV = assign_filters_by_flux(tab1, tab2)
    else:
        # Priorität 2: Command-line Argumente
        label1, label2 = args.filters
        print(f"\nKeine Filter-Info in FITS-Headern.")
        print(f"Verwende --filters Argumente: '{label1}' und '{label2}'")
        
        if label1.upper() == 'B' and label2.upper() == 'V':
            tabB = tab1
            tabV = tab2
            print(f"  → Zuordnung: Bild 1 = B, Bild 2 = V")
        elif label1.upper() == 'V' and label2.upper() == 'B':
            tabB = tab2
            tabV = tab1
            filters_swapped = True
            print(f"  → Zuordnung: Bild 1 = V, Bild 2 = B (vertauscht)")
        else:
            print(f"  Unklare Labels. Verwende Flux-basierte Heuristik...")
            tabB, tabV = assign_filters_by_flux(tab1, tab2)
    
    # WICHTIG: Wenn Filter vertauscht wurden, müssen auch die matched Indizes vertauscht werden!
    # matched enthält (index_in_tab1, index_in_tab2), aber jetzt ist tabB=tab2 und tabV=tab1
    # Also müssen wir zu (index_in_tab2, index_in_tab1) ändern = (index_in_tabB, index_in_tabV)
    if filters_swapped:
        matched = [(iV, iB) for (iB, iV) in matched]
        print(f"  → matched Indizes vertauscht: jetzt (index_in_tabB, index_in_tabV)")
    
    print()  # Leerzeile

    if len(matched) > 0:
        # Zeige Anzahl der Gesamtquellen für Kontext
        n_total_B = len(tab1) if tab1 is not None else 0
        n_total_V = len(tab2) if tab2 is not None else 0

        if args.compare_cmd_methods and not (args.calibrate or args.calibrate_file):
            print("\n=== CMD-Vergleich: aperture vs psf ===")
            comparison_data = {}

            # Bereits berechnete Methode direkt verwenden
            comparison_data[method] = compute_instrumental_cmd_points(
                matched, tabB, tabV, outlier_sigma=args.outlier_sigma
            )

            for compare_method in ('aperture', 'psf'):
                if compare_method in comparison_data:
                    continue

                try:
                    compare_result1 = detect_and_photometry_wrapper(
                        datas[0], method=compare_method, fwhm=args.fwhm,
                        threshold_sigma=args.threshold, aperture_radius=args.aperture
                    ) if datas[0] is not None else (None, None)
                    compare_result2 = detect_and_photometry_wrapper(
                        datas[1], method=compare_method, fwhm=args.fwhm,
                        threshold_sigma=args.threshold, aperture_radius=args.aperture
                    ) if datas[1] is not None else (None, None)
                except Exception as e:
                    print(f"  WARNUNG: {compare_method}-Photometrie fehlgeschlagen: {e}")
                    comparison_data[compare_method] = None
                    continue

                compare_tab1, _ = compare_result1
                compare_tab2, _ = compare_result2
                if compare_tab1 is None or compare_tab2 is None:
                    print(f"  WARNUNG: Keine Daten fuer {compare_method}-Photometrie verfuegbar.")
                    comparison_data[compare_method] = None
                    continue

                compare_matches, _ = match_sources(
                    compare_tab1, compare_tab2, max_sep=args.maxsep,
                    filter_radius=cluster_radius, center=cluster_center,
                    annulus_inner=annulus_inner
                )

                if filters_swapped:
                    compare_matches = [(iV, iB) for (iB, iV) in compare_matches]
                    compare_tabB = compare_tab2
                    compare_tabV = compare_tab1
                else:
                    compare_tabB = compare_tab1
                    compare_tabV = compare_tab2

                comparison_data[compare_method] = compute_instrumental_cmd_points(
                    compare_matches, compare_tabB, compare_tabV,
                    outlier_sigma=args.outlier_sigma
                )

            plot_cmd_method_comparison(
                comparison_data.get('aperture'),
                comparison_data.get('psf'),
                selected,
                save=save_cmd_compare_path,
                calibrated=False,
                dereddened=False
            )
        
        # Kalibrierung mit Referenzstern falls gewünscht
        calibrated_mags_B = None
        calibrated_mags_V = None
        calibrated_colors = None
        
        if args.calibrate or args.calibrate_file:
            print("\n=== Referenzstern-Kalibrierung ===")
            
            # Lade Kalibrierung aus Datei falls angegeben und --calibrate nicht gesetzt
            if args.calibrate_file and not args.calibrate:
                ref_index_B, mag_B_true, mag_V_true, ref_x, ref_y = load_calibration(
                    calibrate_file_path, [tab1, tab2], matched
                )
            else:
                # Interaktive Auswahl
                # Wenn Filter vertauscht wurden, müssen auch die Bilder vertauscht werden!
                images_for_calib = [images[1], images[0]] if filters_swapped else images
                ref_index_B, mag_B_true, mag_V_true, ref_x, ref_y = select_reference_star(
                    images_for_calib, [tabB, tabV], matched, aperture_for_plot, calibrate_file_path, folder
                )
            
            if ref_index_B is not None:
                calibrated_mags_B, calibrated_mags_V, calibrated_colors, matched_calibrated = \
                    calibrate_magnitudes(matched, tabB, tabV, ref_index_B, mag_B_true, mag_V_true)
                
                if calibrated_mags_B is not None:
                    print(f"Kalibrierung erfolgreich! {len(calibrated_mags_B)} Sterne kalibriert.")
                    
                    # Speichere Photometrie-Tabelle
                    # Hole Koordinaten der kalibrierten Sterne (aus tabB, da matched_calibrated auf tabB basiert)
                    x_coords = []
                    y_coords = []
                    for idx_b, idx_v in matched_calibrated:
                        x_coords.append(tabB['x'][idx_b])
                        y_coords.append(tabB['y'][idx_b])
                    
                    phot_table_path = get_output_path(folder, 'photometry.txt')
                    save_photometry_table(phot_table_path, calibrated_mags_V, calibrated_colors,
                                        x_coords=np.array(x_coords), y_coords=np.array(y_coords))
                    
                    # Isochronen-Anpassung falls gewünscht
                    if args.fit_isochrone:
                        # Frage nach Rötung falls nicht angegeben
                        ebv = args.ebv
                        if ebv is None:
                            try:
                                ebv_input = input("\nRötung E(B-V) in mag (z.B. 0.05): ")
                                ebv = float(ebv_input)
                            except (ValueError, EOFError):
                                print("Ungültige Eingabe. Überspringe Isochronen-Anpassung.")
                                ebv = None
                        
                        if ebv is not None:
                            result = fit_isochrone_distance_age(
                                calibrated_colors, calibrated_mags_V, ebv,
                                age_range=tuple(args.age_range),
                                output_folder=folder,
                                fixed_age_myr=args.age
                            )
                            
                            if result is not None:
                                # Entpacke Ergebnisse
                                best_age = result['age_log']
                                dm = result['dm']
                                distance = result['distance_pc']
                                chi2 = result['chi2']
                                
                                # Speichere Ergebnisse
                                results_path = get_output_path(folder, 'cluster_parameters.txt')
                                try:
                                    with open(results_path, 'w') as f:
                                        f.write("# Sternhaufen-Parameter\n")
                                        f.write(f"# Erstellt am: {np.datetime64('now')}\n")
                                        f.write("#\n")
                                        f.write(f"E(B-V) = {ebv:.3f} mag\n")
                                        f.write(f"log(Alter/Jahre) = {best_age:.2f}\n")
                                        f.write(f"Alter = {10**best_age / 1e6:.1f} Millionen Jahre\n")
                                        f.write(f"Entfernungsmodul = {dm:.2f} mag\n")
                                        f.write(f"Entfernung = {distance:.0f} pc = {distance/1000:.2f} kpc\n")
                                        f.write(f"Chi²/N = {chi2:.2f}\n")
                                    print(f"Parameter gespeichert: {results_path}")
                                except Exception as e:
                                    print(f"WARNUNG: Konnte Parameter nicht speichern: {e}", file=sys.stderr)
                                
                                # Bereite Isochrone-Parameter für finales CMD vor
                                # Wenn Isochrone vorhanden, zeige dereddened Daten
                                isochrone_params = result
                                # Dereddene die Daten für das finale CMD
                                calibrated_colors_dered = calibrated_colors - ebv
                                calibrated_mags_V_dered = calibrated_mags_V - 3.1 * ebv
                                calibrated_mags_B_dered = calibrated_mags_B - 4.1 * ebv
                            else:
                                isochrone_params = None
                                calibrated_colors_dered = calibrated_colors
                                calibrated_mags_V_dered = calibrated_mags_V
                                calibrated_mags_B_dered = calibrated_mags_B
                        else:
                            isochrone_params = None
                            calibrated_colors_dered = calibrated_colors
                            calibrated_mags_V_dered = calibrated_mags_V
                            calibrated_mags_B_dered = calibrated_mags_B
                    else:
                        isochrone_params = None
                        calibrated_colors_dered = calibrated_colors
                        calibrated_mags_V_dered = calibrated_mags_V
                        calibrated_mags_B_dered = calibrated_mags_B

                    if args.compare_cmd_methods:
                        print("\n=== CMD-Vergleich: aperture vs psf (kalibriert) ===")
                        comparison_data = {}

                        compare_ebv = None
                        if isochrone_params is not None and 'ebv' in isochrone_params:
                            compare_ebv = isochrone_params['ebv']
                        elif args.ebv is not None:
                            compare_ebv = args.ebv

                        if compare_ebv is None:
                            print("  Hinweis: Kein E(B-V) gesetzt, Vergleich wird kalibriert aber nicht dereddened geplottet.")

                        ref_match_dist = args.maxsep * 3.0
                        if aperture_for_plot is not None:
                            ref_match_dist = max(ref_match_dist, aperture_for_plot * 2.0)

                        for compare_method in ('aperture', 'psf'):
                            try:
                                compare_result1 = detect_and_photometry_wrapper(
                                    datas[0], method=compare_method, fwhm=args.fwhm,
                                    threshold_sigma=args.threshold, aperture_radius=args.aperture
                                ) if datas[0] is not None else (None, None)
                                compare_result2 = detect_and_photometry_wrapper(
                                    datas[1], method=compare_method, fwhm=args.fwhm,
                                    threshold_sigma=args.threshold, aperture_radius=args.aperture
                                ) if datas[1] is not None else (None, None)
                            except Exception as e:
                                print(f"  WARNUNG: {compare_method}-Photometrie fehlgeschlagen: {e}")
                                comparison_data[compare_method] = None
                                continue

                            compare_tab1, _ = compare_result1
                            compare_tab2, _ = compare_result2
                            if compare_tab1 is None or compare_tab2 is None:
                                print(f"  WARNUNG: Keine Daten fuer {compare_method}-Photometrie verfuegbar.")
                                comparison_data[compare_method] = None
                                continue

                            compare_matches, _ = match_sources(
                                compare_tab1, compare_tab2, max_sep=args.maxsep,
                                filter_radius=cluster_radius, center=cluster_center,
                                annulus_inner=annulus_inner
                            )

                            if filters_swapped:
                                compare_matches = [(iV, iB) for (iB, iV) in compare_matches]
                                compare_tabB = compare_tab2
                                compare_tabV = compare_tab1
                            else:
                                compare_tabB = compare_tab1
                                compare_tabV = compare_tab2

                            ref_idx_compare = find_nearest_source_index(
                                compare_tabB, ref_x, ref_y, max_distance=ref_match_dist
                            )
                            if ref_idx_compare is None:
                                print(f"  WARNUNG: Referenzstern fuer {compare_method} nicht gefunden.")
                                comparison_data[compare_method] = None
                                continue

                            comparison_data[compare_method] = compute_calibrated_cmd_points(
                                compare_matches,
                                compare_tabB,
                                compare_tabV,
                                ref_idx_compare,
                                mag_B_true,
                                mag_V_true,
                                outlier_sigma=args.outlier_sigma,
                                ebv=compare_ebv
                            )

                        plot_cmd_method_comparison(
                            comparison_data.get('aperture'),
                            comparison_data.get('psf'),
                            selected,
                            save=save_cmd_compare_path,
                            calibrated=True,
                            dereddened=(compare_ebv is not None)
                        )
                    
                    # Zeige kalibriertes Diagramm (dereddened wenn Isochrone vorhanden)
                    plot_color_magnitude(matched_calibrated, tabB, tabV, selected, 
                                       n_total_B=n_total_B, n_total_V=n_total_V,
                                       calibrated=True, 
                                       mags_B=calibrated_mags_B_dered, 
                                       mags_V=calibrated_mags_V_dered,
                                       colors=calibrated_colors_dered,
                                       save=save_cmd_path,
                                       isochrone_params=isochrone_params,
                                       outlier_sigma=args.outlier_sigma)
                else:
                    print("Kalibrierung fehlgeschlagen. Zeige instrumentale Magnituden.")
                    plot_color_magnitude(matched, tabB, tabV, selected, 
                                       n_total_B=n_total_B, n_total_V=n_total_V, 
                                       save=save_cmd_path,
                                       outlier_sigma=args.outlier_sigma)
            else:
                print("Keine Kalibrierung durchgeführt. Zeige instrumentale Magnituden.")
                plot_color_magnitude(matched, tabB, tabV, selected, 
                                   n_total_B=n_total_B, n_total_V=n_total_V, 
                                   save=save_cmd_path,
                                   outlier_sigma=args.outlier_sigma)
        else:
            # Keine Kalibrierung - zeige instrumentale Magnituden
            plot_color_magnitude(matched, tabB, tabV, selected, 
                               n_total_B=n_total_B, n_total_V=n_total_V, 
                               save=save_cmd_path,
                               outlier_sigma=args.outlier_sigma)
    else:
        print("Keine übereinstimmenden Quellen gefunden für das Color-Magnitude-Diagramm.")
    
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
