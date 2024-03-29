% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/REWarDS.R
\name{definedDailyDose}
\alias{definedDailyDose}
\title{Defined daily dose}
\usage{
definedDailyDose(
  data,
  WHO_ddd,
  dspd_qty,
  strength,
  id,
  tot_dose_disp = NULL,
  Pt_level = FALSE
)
}
\arguments{
\item{data}{Sample simulated data.
Data may have multiple rows per person (one row per prescription fill).
Required columns include:
1. ID: Patient's unique identification number
2. ServDate: Date on which each prescription was filled.
3. DSPD_QTY: Dispensed quantity: Number of tablets dispensed to patient at each prescription fill.
4. strength: Strength of the tablets dispensed.}

\item{WHO_ddd}{The defined daily dose in mg as specified by the WHO.}

\item{dspd_qty}{Dispensed quantity: Number of the dispensed tablets to the patient
at each prescription fill.}

\item{strength}{Strength of the tablet dispensed in milligrams.}

\item{id}{Unique patient identification number.}

\item{tot_dose_disp}{Total dose dispensed:
dispensed quantity x strength of the tablets dispensed for
each prescription fill.}

\item{Pt_level}{When TRUE, the estimated daily dose and days' supply are averaged
for the patient.}
}
\value{
definedDailyDose returns a dataset called "DDD_result". This data set includes all the
         variables originally in the data, plus the following:

DDD_Rx_dose:  Daily dose for prescription.

DDD_Rx_DS:    Days' supply for prescription.

DDD_Pt_dose:  Average daily dose for patient.

DDD_Pt_DS:    Average days' supply for patient.
}
\description{
Uses the World Health Organization's (WHO) defined daily dose method to compute the daily
dose and days' supply for prescriptions. This method assumes an average daily consumption
of a fixed dose, the defined daily dose (DDD), specified by WHO on their website at
https://www.whocc.no/atc_ddd_index/?code=B01AA03.
}
\details{
The DDD method can be used for any medication.
However, its accuracy has been shown to differ between drug classes.
}
\examples{
#Patient collects 100 tablets of 5 mg warfarin on January 3rd,
#and 100 tablets of 7 mg warfarin on February 1st.

#Generate a simulated dataset

library(dplyr)
n_patients <- 10
n_records <- 80
data <- data.frame(ID = rep(c(1 : n_patients), each = n_records))
data \%>\%
  group_by(ID) \%>\%
  mutate(ServDate = as.Date('2020/01/01') + abs(round(rnorm(n = 80, 700, 330))),
         DSPD_QTY = abs(round(rnorm(n = 80, 43, 28))),
         strength = abs(round(rnorm(n = 80, 4, 1))))  -> data
data <- as.data.frame(data)

# Using 1 unit DDD:

data_new <- definedDailyDose (data, WHO_ddd = 7.5, Pt_level = TRUE,
                              id = "ID",dspd_qty = "DSPD_QTY", strength = "strength",
                              tot_dose_disp = NULL)

#WHO_ddd is set as 7.5 mg as that is the defined daily dose set by WHO for warfarin.

#tot_dose_disp: 500mg on January 3rd and 700mg for February 1st.
#DDD_Rx_dose: 7.5 mg for each prescription fill
#DDD_Rx_DS is: For Jan 3rd:  500/7.5 = 66.66 day;
#              For Feb 1st: 700/7.5=93.33 days

#Pt_level can be set as TRUE to get mean values for each patient
#DDD_Pt_dose: (7.5+ 7.5)/2 = 7.5 mg
#DDD_Pt_DS: (66.66+ 93.33)/2 = 79.99 days


}
