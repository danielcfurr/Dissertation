// Find r2 for values of tau ---------------------------------------------------

// Get fixed item effects variance
run sim_dataset.do
sim_dataset, tau(0)
quietly summarize xB
local var_xB = r(Var)

// Make dataset of possible values of tau and r-square
clear
set obs 151
egen double tau = fill(0(.01)1.5)
replace tau = round(tau, .01)
generate tausq = tau^2
generate rsq = `var_xB' / (`var_xB' + tausq)
generate selected = inlist(tau, 0, .1, .5, 1)

// Plot values of r-square against tau
twoway line rsq tau || scatter rsq tau if selected, ///
	plotregion(margin(zero)) yscale(range(0 1) noextend) ylabel(0(.2)1) ///
	msymbol(circle_hollow) msize(1.75) legend(off) ///
	graphregion(fcolor(white)) xtitle({&tau}) ytitle(R{superscript:2})
graph export figs/rsq_vs_tau.pdf, replace

// A program to get list of values from a variable and write this to latex
// macro. Macro option {#1} allows inclusion of "and".
capture program drop writelist
program define writelist
	syntax varlist(max=1) [if] [, file(string) macro(string) fmt(string)]
	if ("`fmt'" == "") local fmt = "%9.2f"
	preserve
		if("`if'" != "") quietly keep `if'
		local newvalue = strofreal(`varlist'[1], "`fmt'")
		local list = "`newvalue'"
		forvalues n = 2/`=_N - 1' {
			local newvalue = strofreal(`varlist'[`n'], "`fmt'")
			local list = "`list', `newvalue'"
		}
		local newvalue = strofreal(`varlist'[`=_N'], "`fmt'")
		local list = "`list', {#1} `newvalue'"
	restore
	local command = "\newcommand{\\`macro'}[1][]{`list'}"
	display subinstr("`list'", "{#1} ", "", 1)
	display "`command'"
	if("`file'" != "") file write `file' "`command'" _n
end

// Write a file of latex macros for generating values
file open myfile using "figs\macros.tex", write replace
	local upsilon2 = strofreal(`var_xB', "%9.2f")
	file write myfile "\newcommand{\genupsilonsq}{`upsilon2'}" _n
	writelist tau if selected, file(myfile) macro(gentau) fmt("%9.2f")
	writelist tausq if selected, file(myfile) macro(gentausq) fmt("%9.2f")
	writelist rsq if selected, file(myfile) macro(rsq) fmt("%9.2f")
file close myfile


// Plot N effective parameters -------------------------------------------------

clear
local sims : dir . files "sim1_*.dta"
foreach f in `"`sims'"' {
   append using `f'
}
drop if insample == . // Not usually necessary

foreach var of varlist aic-newboth {
	generate p_`var'_ = (`var' - insample) / 2
}

collapse (mean) p_*_, by(model tau)

// Write values of AIC and BIC to latex macros.
file open myfile using "figs\macros.tex", write append
	writelist p_aic_ if tau == 0, file(myfile) macro(aic) fmt("%9.2f")
	writelist p_bic_ if tau == 0, file(myfile) macro(bic) fmt("%9.2f")
file close myfile

reshape wide p_*_, i(tau) j(model)

// A program to plot effective n paramters
capture program drop p_plot
program define p_plot
	syntax varlist, [save(string)]
	twoway connect `varlist' tau, ///
		msymbol(circle triangle square) ///
		graphregion(fcolor(white)) ///
		xtitle({&tau}) ytitle("Effective N Parameters") ///
		legend(label(1 "Model 1") label(2 "Model 2") label(3 "Model 3") ///
			rows(1) region(lcolor(white)))	
	if "`save'" != "" {
		graph export "`save'", replace
	}
end

p_plot p_aic_*
p_plot p_bic_*
p_plot p_newpersons_*, save("figs/p_newpersons.pdf")
p_plot p_newitems_*, save("figs/p_newitems.pdf")
p_plot p_newboth_*, save("figs/p_newboth.pdf")


// Plot selection results ------------------------------------------------------

clear
local sims : dir . files "sim1_*.dta"
foreach f in `"`sims'"' {
   append using `f'
}
drop if insample == . // Not usually necessary

reshape wide insample-newboth, i(tau seed) j(model)

// Variables to indicate which model selected for each metric
foreach stub in aic bic newpersons newitems newboth {
	egen double min_`stub' = rowmin(`stub'*)
	forvalues m = 1/3 {
		generate select_`stub'_`m' = (`stub'`m' == min_`stub') if (`stub'`m' < .)
	}
}
drop min_*

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
			graphregion(fcolor(white)) ///
			legend(rows(1) label(1 "Model 1") label(2 "Model 2") ///
				label(3 "Model 3") region(lcolor(white))) ///
			bar(1, fintensity(*1)) bar(2, fintensity(*.7)) bar(3, fintensity(*.4)) ///
			ytitle(Proportion of times selected) 
	restore
	if "`save'" != "" {
		graph export "`save'", replace
	}
end

selectbar select_aic_*, save("figs/select_aic.pdf")
selectbar select_bic_*, save("figs/select_bic.pdf")
selectbar select_newpersons_*, save("figs/select_newpersons.pdf")
selectbar select_newitems_*, save("figs/select_newitems.pdf")
selectbar select_newboth_*, save("figs/select_newboth.pdf")


// Plot N effective parameters for CV over items -------------------------------

clear
local sims : dir . files "sim2_*.dta"
foreach f in `"`sims'"' {
   append using `f'
}
drop if insample == . // Not usually necessary

foreach var of varlist aic-newitems {
	generate p_`var'_ = (`var' - insample) / 2
}

collapse (mean) p_*_, by(model tau)

// Write values of AIC and BIC to latex macros.
file open myfile using "figs\macros.tex", write append
	writelist p_aic_ if tau == 1, file(myfile) macro(aicitem) fmt("%9.2f")
	//writelist p_bic_ if tau == 0, file(myfile) macro(bic) fmt("%9.2f")
file close myfile

reshape wide p_*_, i(tau) j(model)

p_plot p_aic_*
p_plot p_newitems_*, save("figs/p_newitems.pdf")


// Plot selection results for CV over items ------------------------------------

clear
local sims : dir . files "sim2_*.dta"
foreach f in `"`sims'"' {
   append using `f'
}
drop if insample == . // Not usually necessary

reshape wide insample-newitems, i(tau seed) j(model)

// Variables to indicate which model selected for each metric
foreach stub in aic bic newitems {
	egen double min_`stub' = rowmin(`stub'*)
	forvalues m = 1/3 {
		generate select_`stub'_`m' = (`stub'`m' == min_`stub') if (`stub'`m' < .)
	}
}
drop min_*

selectbar select_aic_*, save("figs/select_aic2.pdf")
//selectbar select_bic_*, save("figs/select_bic2.pdf")
selectbar select_newitems_*, save("figs/select_newitems2.pdf")

