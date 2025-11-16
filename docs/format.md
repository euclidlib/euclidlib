# 📁 Format

## Are Euclid products homogeneous?

Original Euclid Consortium Science Ground Segment (SGS) products are not homogeneous beyond being delivered as FITS files. For Level 3 (LE3) products, `euclidlib` performs homogenisation when reading them: it maps heterogeneous FITS inputs to the same internal dataclass formats, enforces the spin/component array shapes and binning conventions described below, and extracts a consistent set of metadata (scale arrays, limits, weights, and header fields). These conversions allow downstream code to treat LE3 products uniformly even if their raw origins differ.

## Photometric products

`euclidlib` reads photometric redshift distributions and returns a canonical representation:

- A 1D redshift array.
- A mapping (Python dictionary) from redshift‑bin index → $$n(z)$$ array

## Level 3 (LE3) Cosmology-Ready Products

These products correspond to measurements of weak lensing (cosmic shear)
and galaxy clustering derived from the Euclid catalogues.

### 1. Spin Structure

Euclid observables map to different array shape. For this, they are grouped by spin; the array shape follows the number of components of each field. Think of each field as a small vector (components ordered consistently) and correlations as outer products of those vectors.

| Spin | Field                    | Components | Component order | Auto-correlation shape |
| ---- | ------------------------ | ---------: | --------------- | ---------------------: |
| 0    | POS (position / density) |          1 | [pos]           |                (1 × 1) |
| 2    | SHE (shear)              |          2 | [e1, e2]        |                (2 × 2) |

Cross-correlations obey the same rule: shape = (components of field A) × (components of field B).

Examples:

```
('SHE', 'SHE', 1, 1)  → shape = (2, 2) # first component is always E-mode
('POS', 'SHE', 1, 2)  → shape = (1, 2)
```

When loaded with `euclidlib`, products are returned as a nested dictionary keyed by:
field1 → field2 → redshift bin 1 → redshift bin 2 → measurement array (with the shapes above). Bins act as independent fields, so every bin pair stores its own component-shaped array.

### 2. `cosmolib` dataclasses

Each measurement is encapsulated in a dedicated Python dataclass corresponding to its observable. For example:

```
AngularPowerSpectra(a, ell=np.(array), axis=scale_axis, ...)
```

These Python dataclasses are provided by the lightweight (`cosmolib` package)[https://github.com/astro-ph/cosmolib]. This format is also used by other python packages such as `heracles`, `cloelib`, `cloelike` and `Spaceborne`.

### 3. Attributes

Each entry of the dictionary contains a dataclass with attributes:

- The array containing the measurements
- The scale ($\theta$ or $\ell$)
- Upper and lower limits of that scale
- Measurement weights
- Additional metadata extracted from the FITS headers

### 4. Real vs Harmonic Space

Real-space and harmonic-space measurements share the same underlying spin structure.
For example:

| Real Space                  | Harmonic Space                      |
| --------------------------- | ----------------------------------- |
| $w(\theta)$                 | $C_\ell^{gg}$                       |
| $(\gamma_t, \gamma_\times)$ | $(C_\ell^{gE}, C_\ell^{gB})$        |
| $(\xi_{+}, \xi_{-})$        | $(C_\ell^{EE}, C_\ell^{BB}, \dots)$ |

Because of this shared structure, all these quantities are stored in a uniform way using the `euclidlib` format—even though their original data products are not homogeneous. In practice, `euclidlib` handles the required conversions automatically.
