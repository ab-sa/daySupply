% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/REWarDS.R
\name{fixedTablet}
\alias{fixedTablet}
\title{Fixed tablet}
\usage{
fixedTablet(
  data,
  tablet = 1,
  dspd_qty,
  strength,
  id,
  serv_date,
  tot_dose_disp = NULL,
  Pt_level = FALSE
)
}
\arguments{
\item{data}{Sample simulated data}

\item{tablet}{Number of tablets assumed to be consumed by the patient per day.
Default=1}

\item{dspd_qty}{Dispensed quantity: Number of the dispensed tablets to the patient
at each prescription fill}

\item{strength}{Strength of the tablet dispensed in milligrams}

\item{id}{Unique patient identification number}

\item{tot_dose_disp}{Total dose dispensed:
dispensed quantity x strength of the tablets dispensed for
each prescription fill}

\item{Pt_level}{When TRUE, the estimated daily dose and days' supply are averaged
for the patient}

\item{servDate}{Date of the prescription fill}
}
\value{
tot_dose_disp:            Total dose dispensed at prescription fill:
                                  dispensed quantity x strength of the tablet dispensed

fixed_1_tab_Rx_dose:      Daily dose for prescription

fixed_1_tab_Rx_DS:        Days' supply for prescription

fixed_1_tab_Pt_dose:      Average daily dose for patient

fixed_1_tab_Pt_DS:        Average days' supply for patient
}
\description{
fixedTablet is a function that computes the daily dose and days'
supply for prescriptions by assuming an average daily consumption of a
fixed number of tablets (usually 1) per day by the patient.
}
\details{
The fixed tablet method can be used for any medication.
However, its accuracy has been shown to differ between drug classes.
}
\examples{
#Patient collects 100 tablets of 5 mg warfarin  on January 3rd,
#and 100 tablets of 7 mg warfarin on February 1st.

#Assuming consumption of 1 tablet per day:
fixedTablet(data = data, tablet = 1, Pt_level = F, id = "ID",
            dspd_qty = "DSPD_QTY", strength = "strength",
            serv_date = "ServDate", tot_dose_disp = NULL)

#tot_dose_disp= 500mg on January 3rd and
                700 mg for February 1st.
#fixed_1_tab_Rx_dose: 5 mg for the prescription refill on Jan 3rd, 7 mg for prescription refill on Feb 1st.
#fixed_1_tab_Rx_DS is: For Jan 3rd:  500/5= 100 day;  For Feb 1st: 700/7= 100 days

#pt_level can be set as TRUE to get mean values for each patient
#DDD_1_Pt_dose: (5+ 7)/2 = 6 mg
#DDD_1_Pt_DS: (100+100)/2 = 100 days



}