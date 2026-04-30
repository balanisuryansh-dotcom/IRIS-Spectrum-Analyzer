import numpy as np
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output
from astropy.io import fits
from scipy.optimize import curve_fit
import sys

# Constants
SPEED_OF_LIGHT = 299792.458  # km/s


def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-(x - mean) ** 2 / (2 * stddev ** 2))


def run_interactive_spectrum_viewer(fits_path):
    try:
        print(f"Loading FITS file from: {fits_path}")

        with fits.open(fits_path) as hdul:
            if len(hdul) < 2 or hdul[1].data is None:
                raise ValueError("Invalid FITS structure: missing data in HDU[1]")

            cube = hdul[1].data
            header = hdul[1].header

            crval1 = header.get('CRVAL1')
            cdelt1 = header.get('CDELT1')
            crpix1 = header.get('CRPIX1', 1.0)

            if crval1 is None or cdelt1 is None:
                raise ValueError("Missing wavelength calibration (CRVAL1/CDELT1) in FITS header.")

        num_wavelengths = cube.shape[0]
        wavelength = crval1 + cdelt1 * (np.arange(num_wavelengths) - (crpix1 - 1))

        print(f"Wavelength range: {wavelength.min():.2f} – {wavelength.max():.2f} Å")

        mean_intensity = np.mean(cube, axis=0)

        app = Dash(__name__)

        app.layout = html.Div([
            html.H2("Interactive IRIS Spectrum Viewer"),

            dcc.Graph(
                id='image-plot',
                figure=go.Figure(
                    data=go.Heatmap(
                        z=mean_intensity,
                        colorscale='Viridis',
                        colorbar=dict(title='Avg Intensity'),
                        hovertemplate='X: %{x}<br>Y: %{y}<extra></extra>'
                    ),
                    layout=go.Layout(
                        title="Click a Pixel to View Spectrum",
                        xaxis_title="X",
                        yaxis_title="Y"
                    )
                )
            ),

            html.Div([
                html.Label("X:"), dcc.Input(id='x-input', type='number', value=0, min=0, max=cube.shape[2] - 1),
                html.Label("Y:"), dcc.Input(id='y-input', type='number', value=0, min=0, max=cube.shape[1] - 1),
            ], style={'margin': '10px', 'display': 'flex', 'gap': '10px'}),

            html.Div([
                html.Label("Min WL (Å):"),
                dcc.Input(id='min-wl', type='number', value=float(np.min(wavelength))),
                html.Label("Max WL (Å):"),
                dcc.Input(id='max-wl', type='number', value=float(np.max(wavelength))),
            ], style={'margin': '10px', 'display': 'flex', 'gap': '10px'}),

            html.Div([
                html.Label("Rest Wavelength Mode:"),
                dcc.Dropdown(
                    id='rest-wl-mode',
                    options=[
                        {'label': 'Default', 'value': 'default'},
                        {'label': 'Si IV (1393 Å)', 'value': 'siiv_1393'},
                        {'label': 'Si IV (1402 Å)', 'value': 'siiv_1402'},
                        {'label': 'C II (1335 Å)', 'value': 'cii_1335'},
                        {'label': 'Mg II k (2796 Å)', 'value': 'mgii_2796'}
                    ],
                    value='default',
                    clearable=False
                )
            ], style={'margin': '10px', 'width': '50%'}),

            html.Div([
                html.Label("Show Fits:"),
                dcc.Checklist(
                    id='fit-toggle',
                    options=[
                        {'label': 'Curve Fit', 'value': 'curve'},
                        {'label': 'FWHM Fit', 'value': 'fwhm'}
                    ],
                    value=[]
                ),
            ], style={'margin': '10px'}),

            dcc.Graph(id='spectrum-plot')
        ])

        @app.callback(
            Output('spectrum-plot', 'figure'),
            [Input('x-input', 'value'),
             Input('y-input', 'value'),
             Input('min-wl', 'value'),
             Input('max-wl', 'value'),
             Input('rest-wl-mode', 'value'),
             Input('fit-toggle', 'value')]
        )
        def update_spectrum(x_input, y_input, min_wl, max_wl, rest_wl_mode, toggles):
            try:
                x, y = int(x_input), int(y_input)

                if not (0 <= y < cube.shape[1] and 0 <= x < cube.shape[2]):
                    return go.Figure().update_layout(title="Coordinates out of bounds")

                spectrum = cube[:, y, x]
                mask = (wavelength >= min_wl) & (wavelength <= max_wl)

                wl = wavelength[mask]
                spec = spectrum[mask]

                fig = go.Figure()
                fig.add_trace(go.Scatter(x=wl, y=spec, mode='lines', name="Spectrum"))

                centroid, sigma_cf, sigma_fwhm = None, None, None

                if 'curve' in toggles:
                    try:
                        p0 = [np.max(spec), wl[np.argmax(spec)], np.std(wl)]
                        popt, _ = curve_fit(gaussian, wl, spec, p0=p0, maxfev=10000)
                        _, centroid, sigma_cf = popt
                        fig.add_trace(go.Scatter(x=wl, y=gaussian(wl, *popt),
                                                 mode='lines', name="Curve Fit"))
                    except Exception as e:
                        print(f"Curve fit failed: {e}")

                if 'fwhm' in toggles:
                    peak_idx = np.argmax(spec)
                    peak_val = spec[peak_idx]
                    half_max = peak_val / 2

                    left, right = peak_idx, peak_idx
                    while left > 0 and spec[left] > half_max:
                        left -= 1
                    while right < len(spec) - 1 and spec[right] > half_max:
                        right += 1

                    if left != right:
                        fwhm = wl[right] - wl[left]
                        sigma_fwhm = fwhm / 2.355

                    centroid = wl[peak_idx] if centroid is None else centroid
                    fig.add_trace(go.Scatter(x=wl,
                                             y=gaussian(wl, peak_val, centroid, sigma_fwhm or 0),
                                             mode='lines',
                                             name="FWHM Fit"))

                if centroid is not None:
                    rest_map = {
                        'default': np.floor(centroid),
                        'siiv_1393': 1393.76,
                        'siiv_1402': 1402.77,
                        'cii_1335': 1335.71,
                        'mgii_2796': 2796.35
                    }

                    rest_wl = rest_map.get(rest_wl_mode, np.floor(centroid))
                    doppler_velocity = ((centroid - rest_wl) / rest_wl) * SPEED_OF_LIGHT

                    fig.update_layout(
                        title=f"Centroid: {centroid:.3f} Å | Velocity: {doppler_velocity:.2f} km/s"
                    )

                fig.update_layout(xaxis_title="Wavelength (Å)", yaxis_title="Intensity")
                return fig

            except Exception as e:
                return go.Figure().update_layout(title=f"Error: {str(e)}")

        print("Starting Dash server...")
        app.run(debug=True, port=8070)

    except Exception as e:
        print(f"Error loading FITS file: {e}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python app.py <path_to_fits>")
    else:
        run_interactive_spectrum_viewer(sys.argv[1])
