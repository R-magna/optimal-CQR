# Optimal Conformalized Quantile Regression

# Abstract:
Uncertainty quantification (UQ) has become a central topic in modern machine learning, as reliable predictive intervals are crucial for decision-making under uncertainty. Conformal Prediction (CP) provides a distribution-free framework for constructing prediction sets with finite-sample coverage guarantees, without any assumptions about the underlying model. However, standard CP methods such as Conformalized Quantile Regression (CQR) typically rely on fixed quantile levels, which can lead to inefficiencies under heteroskedasticity or asymmetric conditional distributions.
In this work, we introduce Optimal Conformalized Quantile Regression (Optimal-CQR), a refinement of CQR that adaptively selects quantile levels based on covariates to minimize interval width while preserving coverage guarantees. We establish its theoretical foundations, proving that Optimal-CQR maintains finite-sample validity and converges to oracle intervals under mild regularity conditions, achieving asymptotic length efficiency compared to standard CQR. Extensive experiments on twelve benchmark datasets and a real-world application to superconductivity data demonstrate that Optimal-CQR can improve efficiency, often reducing interval length by up to 35% without compromising coverage. These results highlight the practical potential of Optimal-CQR for uncertainty-aware predictive modeling in complex domains.

The code and the article are inspired by the thesis 
# "Conformal Prediction for Regression: a proposal to improve Conformalized Quantile Regression"
Candidate: Roberto Magnacavallo
Tutor: Matteo Borrotti
Co-Tutor: Valentina Zangirolami

# Note
This repository contains code and data to replicate the applications. 
