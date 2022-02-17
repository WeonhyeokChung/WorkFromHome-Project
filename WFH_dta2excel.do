********************************************************************************
*
*	Project Name: WFH dta to excel
*	Date: Feb 17th 2022
*	
********************************************************************************

*------------------------------------------------------------------------------*
* 							Prep											   *
*------------------------------------------------------------------------------*

cd "/Users/BrianChung/Dropbox/Github"

use WorkFromHomeProject/treatment_effect, clear
export delimited using WorkFromHomeProject/treatment_effect.csv, replace

	
