# Tesi

## 1 Introduction

This git is the repository for my bachelor degree thesis, which consists of studying a truncated state entering a Mach-Zehnder interferometer, that interferes with different states in different experiments.

In particular, the truncated state is a truncated superposition of Fock states with proper normalized coefficients $\psi_k$.

$$
\ket{\psi} = \sum_{k=0}^N{\psi_k \ket{k}}
$$

During the thesis work I studied this state and its interference with different state in different systems, in particular, the truncated state enters the Mach-Zehnder interferometer along the horizontal mode $\hat{a}$, while the other state enters the system along the vertical mode $\hat{b}$, the states studied are

1. the vacuum state $(\ket{0})$ ;
2. a coherent seed $(\ket{\alpha})$, with $\alpha \in \mathbb{C}$ ;
3. a thermal state, which represents the thermal background $(\hat{\varrho}_{th})$ .

The study consists of estimating the quantum Fisher information for a phase $\phi$ introduced on the state by the evolution through a sample placed along the horizontal mode of the interferometer, in particular, the calculation shows that the quantum Fisher information $H$ does not depend on $\phi$, and can be seen as the variance of the operator $\hat{G}$, defined as follows, computed on the initial global states $\ket{\Psi}$ for cases 1. and 2., while assumes a more complicated shape for the thermal state.

$$
\hat{G} = \frac{1}{2} (\hat{N}_a + \hat{N}_b - \hat{a}^\dagger \hat{b} - \hat{a} \hat{b}^\dagger)
$$

## 2 Truncated state & Vacuum state


