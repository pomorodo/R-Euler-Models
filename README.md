# Euler Approximation simulation

Simulation of allele frequency dynamics under the **k-allele variants diffusion model**, implemented in R as part of a Master's thesis.

---

## Model

The model describes the evolution of allele frequencies $z = (z_1, \dots, z_{k-1})$ on the simplex $\Delta^{k-1}$ via a stochastic differential equation (SDE):

$$dz_i = \sum_j L_{ij}(z) dW_j$$

where the diffusion matrix is:

$$A_{ij}(z) = \frac{1}{2N_e} z_i (\delta_{ij} - z_j), \quad L \text{ is such that } L L^\top = A$$

and $N_e$ is the effective population size.

The $k$-th allele frequency is recovered as $z_k = 1 - \sum_{i=1}^{k-1} z_i$.

---


## Project Structure

```
\textit{work in progress…}
```

---

## Usage

```
\textit{work in progress…}
```

---


## Author

**Aymen** — Master's thesis

---

## License

This repository is private during thesis development. All rights reserved.
