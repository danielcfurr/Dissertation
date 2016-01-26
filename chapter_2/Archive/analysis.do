run sim_dataset.do


// Tables for design matrices --------------------------------------------------

sim_dataset, i_factor(1) p_factor(1)
drop if new_person | new_item

// requires outtable
// ssc install outtable
// outtable is not versatile, so must fixed up "raw" latex output

sort x1-x4
mkmat x1-x4 if person == 1, matrix(X)
matrix list X
outtable using "figs/table_x_raw", mat(X) replace nobox asis  ///
	caption("Items design matrix") norowlab clabel("tab:X")

sort w1-w2
mkmat w1-w2 if item == 1, matrix(W)
matrix list W
outtable using "figs/table_w_raw", mat(W) replace nobox asis  ///
	caption("Persons design matrix") norowlab clabel("tab:W")
	
	
// Find r2 for values of tau ---------------------------------------------------

// Make matrix of values of tau and resulting R^2
matrix R = J(100,2,.)
matrix colnames R = tau r2
local i = 0
forvalues tau = .00(.01)1 {
	local ++i
	sim_dataset, i_factor(1) p_factor(0) tau(`tau')
	if "`v_xB'" == "" {
		quietly summarize xB
		local v_xB = r(Var) * (_N - 1) / _N // Correct sample variance to population
	}
	local r2 = `v_xB' / (`tau'^2 + `v_xB')
	matrix R[`i',1] = `tau'
	matrix R[`i',2] = `r2'
}

// Plot R^2 against tau
clear
svmat R, names(col)
generate selected = inlist(tau, float(.0), float(.1), float(.3), float(.5))
twoway line r2 tau || scatter r2 tau if selected, ///
	msymbol(circle_hollow) msize(1.75) legend(off) ///
	graphregion(fcolor(white)) xtitle({&tau}) ytitle(R{superscript:2})
graph export figs/rsq_vs_tau.pdf, replace

// Make a file for latex macros related to selected values of tau and R^2
mkmat tau r2 if selected, matrix(S)
local R = rowsof(S)
forvalues r = 1/`R' {
	display `r'
	local tau_r = strofreal(el(S, `r', 1), "%9.2f")
	local tau2_r = strofreal(el(S, `r', 1)^2, "%9.2f")
	local r2_r = strofreal(el(S, `r', 2), "%9.2f")
	if `r' == 1 {
		local tau_vec = "\{`tau_r'"
		local tau2_vec = "\{`tau2_r'"
		local r2_vec = "\{`r2_r'"
	}
	else if `r' == `R' {
		local tau_vec = "`tau_vec', `tau_r'\}"
		local tau2_vec = "`tau2_vec', `tau2_r'\}"
		local r2_vec = "`r2_vec', `r2_r'\}"
	}
	else {
		local tau_vec = "`tau_vec', `tau_r'"
		local tau2_vec = "`tau2_vec', `tau2_r'"
		local r2_vec = "`r2_vec', `r2_r'"
	}
}
file open myfile using "figs\macros.tex", write replace
local upsilon2 = strofreal(`v_xB', "%9.2f")
file write myfile "\newcommand{\genupsilonsq}{`upsilon2'}" _n
file write myfile "\newcommand{\gentau}{`tau_vec'}" _n
file write myfile "\newcommand{\gentausq}{`tau2_vec'}" _n
file write myfile "\newcommand{\genrsq}{`r2_vec'}" _n
file close myfile


// LLTM ------------------------------------------------------------------------

clear
local sims : dir . files "sim_*.dta"
foreach f in `"`sims'"' {
   append using `f'
}

// Check for convergence of all simulations
tabstat *_converged_*, stat(count sum)
drop *_converged_*

generate n =_n
generate P = "P = " + string(pfactor * 4, "%04.0f")
generate I = "I = " + string(ifactor * 8, "%03.0f")

// Variables to indicate which model selected for each metric
foreach stub in p_dev_01 p_dev_11 p_aic i_dev_11 i_aic {
	egen double min_`stub' = rowmin(`stub'*)
	forvalues m = 1/3 {
		generate select_`stub'_`m' = (`stub'_`m' == min_`stub') if (`stub'_`m' < .)
	}
}
drop min_*

// Variable to indicate which model selected based on p-values
forvalues p = 2/3 {
	generate temp`p' = p_pcheck_`p' < .05
}
generate select_pcheck_1 = !temp2
generate select_pcheck_2 = temp2 & !temp3
generate select_pcheck_3 = temp2 & temp3
drop p_pcheck_1 temp*

// Program to make bar graphs of selected models
capture program drop selectbar
program define selectbar
	syntax varlist [, save(string)]
	tempvar nmissing
	egen `nmissing' = rowmiss(`varlist')
	quietly levelsof tau if `nmissing' == 0, local(taus)
	preserve
		keep if `nmissing' == 0
		quietly levelsof tau, local(taus)
		local ntaus : word count `taus'
		forvalues t = 1/`ntaus' {
			local val : word `t' of `taus'
			local relabel = `"`relabel' "' + `"`t' "{&tau} = `val'""' 
		}
		graph bar `varlist' if `nmissing' == 0, nofill ///	
			over(tau, label(labsize(small)) relabel(`relabel')) ///
			by(P I, note("") graphregion(fcolor(white)) rows(2)) ///
			legend(rows(1) label(1 "Model 1") label(2 "Model 2") ///
				label(3 "Model 3") region(lcolor(white))) ///
			bar(1, fintensity(*1)) bar(2, fintensity(*.7)) bar(3, fintensity(*.4)) ///
			ytitle(Proportion of times selected) 
	restore
	if "`save'" != "" {
		graph export "`save'", replace
	}
end

selectbar select_pcheck*, save("figs/p_pcheck.pdf")
selectbar select_p_aic*, save("figs/p_aic.pdf")
selectbar select_p_dev_01*, save("figs/p_new_person_same_item.pdf")
selectbar select_p_dev_11*, save("figs/p_new_person_new_item.pdf")
selectbar select_i_dev_11*, save("figs/i_new_person_new_item.pdf")
selectbar select_i_aic*, save("figs/i_new_person_new_item.pdf")

preserve
	generate s1 = p_dev_11_1 < p_dev_11_2
	generate s3 = p_dev_11_3 < p_dev_11_2
	collapse (mean) s1 s3, by(tau I P)
	twoway connected s1 s3 tau, yscale(range(0 1) noextend) ///
		by(I P, note("") graphregion(fcolor(white)) rows(2)) ///
		lcolor(navy forest_green) mcolor(navy forest_green) ///
		msymbol(circle triangle) xlabels(0(.1).5) ///
		xtitle({&tau}) ytitle(Proportion of times selected over Model 2) ///
		legend(rows(1) label(1 "Model 1") label(2 "Model 3") region(lcolor(white)))	
	graph export figs/both_line.pdf, replace
restore

