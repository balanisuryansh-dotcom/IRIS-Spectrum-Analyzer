# IRIS Spectrum Analyzer

An interactive web-based tool for visualizing and analyzing spectral data from FITS files.  
Built using Dash and Plotly, this application enables pixel-level spectral inspection, Gaussian fitting, and Doppler velocity estimation in real time.

---

##  Features

- **Pixel-Level Exploration**  
  Select any pixel from the intensity map to view its spectral profile.

- **Interactive Heatmap**  
  Visual representation of average intensity across the dataset.

- **Wavelength Filtering**  
  Focus on specific wavelength ranges for detailed analysis.

- **Gaussian Curve Fitting**  
  Fit spectral lines using a Gaussian model to estimate centroid and spread.

- **FWHM-Based Approximation**  
  Alternative estimation of spectral width using Full Width at Half Maximum.

- **Doppler Velocity Calculation**  
  Compute velocity shifts based on selected or predefined rest wavelengths.

---

## Technologies Used

- Python  
- NumPy  
- Dash  
- Plotly  
- Astropy  
- SciPy  

---

## Input Data

This tool works with spectral data stored in FITS (Flexible Image Transport System) format.

It assumes:
- A 3D data cube (wavelength × Y × X)
- Valid header calibration using:
  - `CRVAL1`
  - `CDELT1`
  - `CRPIX1`

---

## How to Run

### 1. Install dependencies
```bash
pip install numpy dash plotly astropy scipy
