# Conformal_Imputation
This repository contains the script *glucodensities_missing.R* with functions for applying the methodology presented in the manuscript **Personalized Imputation in metric spaces via conformal prediction: Applications in Predicting Diabetes Development with Continuous Glucose Monitoring Information**. It also includes other three files:
* *conformal_imputations_glucodensities.html*, where it is explained step by step, with codes, the application to real data present in this manuscript,
* *conformal_imputations_simulations.html*, where a small simulation is presented,
* *functions*, with some auxiliar functions used in the analyses included in the html files. 

## Personalized Imputation in metric spaces via conformal prediction: Applications in Predicting Diabetes Development with Continuous Glucose Monitoring Information
Marcos Matabuena$^{1,2\*}$, Carla Díaz-Louzao^3, Rahul Gosal^4, Francisco Gude^{5,6,7}

^1 Health Institute of Santiago de Compostela, Spain

^2 Department of Biostatistics, Harvard University, Boston, MA 02115, USA

^3 Department of Mathematics, University of A Coruña, Spain

^4 Department of Epidemiology and Biostatistics, University of South Carolina, USA

^5 ISCIII Support Platforms for Clinical Research, Health Institute of Santiago de Compostela (IDIS), Spain

^6 Concepción Arenal Primary Care Center, Santiago de Compostela, Spain

^7 Department of Psychiatry, Radiology, Public Health, Nursing and Medicine, University of Santiago de Compostela, Spain


\* Corresponding author: mmatabuena@hsph.harvard.edu

The challenge of handling missing data is widespread in modern data analysis, particularly during the preprocessing phase and in various inferential modeling tasks. Although numerous algorithms exist for imputing missing data, the assessment of imputation quality at the patient level often lacks personalized statistical approaches. Moreover, there is a scarcity of imputation methods for metric space based statistical objects. The aim of this paper is to introduce a novel two-step framework that comprises: (i) a imputation methods for statistical objects taking values in metrics spaces, and (ii) a criterion for personalizing imputation using conformal inference techniques. This work is motivated by the need to impute distributional functional representations of continuous glucose monitoring (CGM) data within the context of a longitudinal study on diabetes, where a significant fraction of patients do not have available CGM profiles. The importance of these methods is illustrated by evaluating the effectiveness of CGM data as new digital biomarkers to predict the time to diabetes onset in healthy populations. To address these scientific challenges, we propose: (i) a new regression algorithm for missing responses; (ii) novel conformal prediction algorithms tailored for metric spaces with a focus on density responses within the 2-Wasserstein geometry; (iii) a broadly applicable personalized imputation method criterion, designed to enhance both of the aforementioned strategies, yet valid across any statistical model and data structure. Our findings reveal that incorporating CGM data into diabetes time-to-event analysis, augmented with a novel personalization phase of imputation, significantly enhances predictive accuracy by over ten percent compared to traditional predictive models for time to diabetes.
