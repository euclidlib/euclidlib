# 📁 Data Format

## Are Euclid products homogeneous?

Original Euclid Consortium Science Ground Segment (SGS) products and
other Science Working Group (SWG) products are **not homogeneous**,
beyond the fact that they are delivered as FITS files.

Level 3 (LE3) and other cosmology-ready products (summary statistics)
differ in:

- Internal FITS structure
- Naming conventions
- Array shapes
- Metadata formatting

To ensure consistency as much as possible, `euclidlib` performs **automatic homogenisation
at read time**. It:

- Maps heterogeneous FITS inputs to unified internal dataclass
  formats
- Enforces consistent spin/component array shapes
- Standardises binning conventions
- Extracts a common set of metadata (scales, limits, weights, headers)

This guarantees that downstream cosmological analyses can treat all
products uniformly, regardless of their raw origin, which is highly-useful in likelihood codes.

---

## SGS Photometric Products

`euclidlib` reads photometric redshift distributions and returns a
canonical representation consisting of:

- A 1D redshift array
- A Python dictionary mapping `redshift_bin_index` to a $n(z)$ array

This ensures all photo-z inputs share a consistent internal format.

---

## Level 3 (LE3) and Other Cosmology-Ready Products

### 1. `euclidlib` Output Structure

When a summary-statistics product is loaded, `euclidlib` returns a
**Python dictionary** whose key are tuples:

```python
{
    ('FIELD1', 'FIELD2', zbin1, zbin2): dataclass_instance
}
```

Each key uniquely identifies, for instance:

    (field_1, field_2, redshift_bin_1, redshift_bin_2)

Bins are treated as independent fields, so every bin pair stores its own
component-shaped measurement array.

---

### 2. `cosmolib` Dataclasses

Each measurement is stored in a dedicated Python dataclass corresponding
to its observable.

Example:

```python
AngularPowerSpectra(
    array=np.array,
    ell=np.array,
    scale_axis="ell",
    ...
)
```

These Python dataclasses are provided by the lightweight [cosmolib package](https://github.com/astro-ph/cosmolib). This format is also used by other python packages such as heracles, cloelib, cloelike and Spaceborne.

---

### 3. Spin Structure

Euclid observables are grouped by **spin**, which determines the number
of components and therefore the array shape.

- **POS** = Photometric galaxy positions
- **SHE** = Weak lensing shear
- **SPE** = Spectroscopic galaxy clustering

| Spin | Field                    | Components | Component order | Auto-correlation shape | Quantity stored in attribute |
| ---- | ------------------------ | ---------: | --------------- | ---------------------: | ---------------------------- |
| 0    | POS (position / density) |          1 | [pos]           |                (1 × 1) | array                        |
| 2    | SHE (shear)              |          2 | [e1, e2]        |                (2 × 2) | array                        |
| 0    | SPE (position / density) |          1 | [spe]           |                (1 × 1) | multipoles                   |

Cross-correlations obey the same rule: shape = (components of field A) × (components of field B). Examples:

```python
('SHE', 'SHE', 1, 1).array       → shape = (2, 2) # first component is always E-mode, second is B-mode
# ('SHE', 'SHE', 1, 1).array[0, 0] will have lenght equal to ('SHE', 'SHE', 1, 1).ell
('POS', 'SHE', 1, 2).array       → shape = (1, 2)
# ('SHE', 'SHE', 1, 1).array[0, 0] will have lenght equal to ('SHE', 'SHE', 1, 1).ell
('SPE', 'SPE', 1, 1).multipoles  → shape = (5, Nk)
# ('SPE', 'SPE', 1, 1).multipoles will have lenght
```

---

### 4. Dataclass Attributes

Each dictionary entry contains a dataclass with:

- The measurement array (`array` or `multipoles`)
- The scale (`theta`, `ell`, or `k`)
- Upper and lower scale limits
- Measurement weights
- Metadata extracted from FITS headers

---

### 5. Photometric: Real vs Harmonic Space

Photometric real-space and harmonic-space observables share the same underlying spin
structure.

| Real Space                  | Harmonic Space                      |
| --------------------------- | ----------------------------------- |
| $w(\theta)$                 | $C_\ell^{gg}$                       |
| $(\gamma_t, \gamma_\times)$ | $(C_\ell^{gE}, C_\ell^{gB})$        |
| $(\xi_{+}, \xi_{-})$        | $(C_\ell^{EE}, C_\ell^{BB}, \dots)$ |

Because of this shared structure, all observables are stored in a
uniform internal format by `euclidlib` regardless of the working-space.
