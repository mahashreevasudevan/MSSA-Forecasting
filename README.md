## Multivariate Singular Spectrum Analysis (MSSA) for Dew Point and Relative Humidity short-term forecasting ##

This project uses Multivariate Singular Spectrum Analysis (MSSA) to forecast dew point temperature and relative humidity using historical data. 
It uses the eigendecomposition of a covariance matrix and a linear recurrent forecasting model to predict future values.

## Key Objectives:

- To analyze and forecast short-term, dew point and relative humidity using MSSA.
- To develop a data-driven forecasting model that captures nonlinearities and dependencies between multiple time series.
- To compare reconstructed and original time series for quality assessment of MSSA decomposition.

## Methodology:

1. **Data Acquisition**: Dew point and relative humidity values are read from Excel sheets.
2. **Normalization**: Both time series are normalized (zero mean and unit variance).
3. **MSSA**:
   - Time-delay embedding of each time series.
   - Formation of a joint trajectory matrix.
   - Covariance matrix construction and eigendecomposition.
   - Projection onto principal components (PCs).
4. **Reconstruction**: Selected PCs are used to reconstruct the original time series.
5. **Forecasting**: Future points are predicted using a Linear Recurrent Formula (LRF) derived from dominant eigenvectors.
6. **Denormalization**: Final forecasted values are transformed back to their original scale.

## Model Pipeline:
Excel Input (Dew Point, RH) -> Normalization (Z-score) -> Time-Delayed Embedding -> Covariance Matrix → Eigendecomposition -> Principal Components (PCs) -> Reconstructed Components (RCs) -> Linear Recurrent Forecasting -> Denormalization -> Forecasted Dew Point & Relative Humidity

## Challenges Addressed:

- Modeling nonlinear dependencies between atmospheric variables.
- Noise reduction and component selection for easier reconstructions.
- Building a short-term forecasting model with no prior assumptions on data distribution.

## Results:

- Forecasted both dew point and relative humidity for short-term with low error values.
- Reconstructed time series using dominant PCs.
- Visualized principal components and eigenstructure.
- Exported key matrices, components, and results to Excel.

## Technology and Tools:

- **Language**: MATLAB
- **Libraries/Functions**:
  - `xlsread`, `xlswrite` (Excel I/O)
  - `plot`, `imagesc`, `eig`, `mean`, `std`
- **Techniques**: MSSA, LRF, PCA

## Note:
The dataset used in this project is confidential. Hence, the dataset, the output and plots are not shown.



