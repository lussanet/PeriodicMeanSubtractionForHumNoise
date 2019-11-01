# PeriodicMeanSubtractionForHumNoise
This subtraction filter analyses the shape of the power line interference ("Hum noise") and subtracts it from the time series. Method: moving estimate of Hum noise by periodic median of the high-pass filtered signal. Highly robust to any harmonic contributions and transients in the signal.
