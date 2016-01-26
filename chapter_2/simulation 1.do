run sim_dataset.do
sim_dataset, tau(0)
export delimited using "example_sim_data.csv", replace


// Conduct simulation ----------------------------------------------------------

run sim_dataset.do
local s = 1

import delimited "random.txt", delimiter(tab) varnames(nonames) clear
rename v`s' seed
drop v*
keep in 1/700

egen condition = seq(), from(1) to(7)
generate double tau = .3
replace tau = 0 if condition == 1
replace tau = .1 if condition == 2
replace tau = .3 if condition == 3
replace tau = .5 if condition == 4
generate nitems = 32
replace nitems = 64 if condition == 5
replace nitems = 96 if condition == 6
replace nitems = 128 if condition == 7
generate npersons = 500

generate insample = .
generate aic = .
generate bic = .
generate newpersons = .
generate newitems = .
generate newboth = .
generate failures = 0

mkmat *, matrix(S1)
matrix define S2 = S1
matrix define S3 = S1

local m1 "w1 w2 x2-x4"
local m2 "w1 w2 x2-x4 c.x2#c.x3"
local m3 "w1 w2 x2-x4 c.x2#c.x3 c.x3#c.x4"
local N = _N

timer clear
timer on 1

forvalues i = 1/`N' {

	display "$S_TIME: Started `i' of `N'"

	sim_dataset, seed(`=el(S1, `i', 1)') tau(`=el(S1, `i', 3)') ///
		nitems(`=el(S1, `i', 4)') npersons(`=el(S1, `i', 5)') 

	forvalues m = 1/3 {
	
		local failures = 0
	
		quietly melogit y_newpersons `m`m'' || person:
		if (e(converged) == 0) local ++failures
		matrix B_newpersons = e(b)
		local cols : colfullnames B_newpersons 
		matrix colnames B_newpersons = `=subinstr("`cols'", "_newpersons", "", .)'

		quietly melogit y_newitems `m`m'' || person:
		if (e(converged) == 0) local ++failures
		matrix B_newitems = e(b)
		local cols : colfullnames B_newitems 
		matrix colnames B_newitems = `=subinstr("`cols'", "_newitems", "", .)'
		
		quietly melogit y_newboth `m`m'' || person:
		if (e(converged) == 0) local ++failures
		matrix B_newboth = e(b)
		local cols : colfullnames B_newboth 
		matrix colnames B_newboth = `=subinstr("`cols'", "_newboth", "", .)'

		quietly melogit y `m`m'' || person:
		if (e(converged) == 0) local ++failures
		matrix S`m'[`i', colnumb(S`m', "insample")] = -2 * e(ll)
		
		quietly estat ic
		matrix S`m'[`i', colnumb(S`m', "aic")] = el(r(S), 1, 5)
		matrix S`m'[`i', colnumb(S`m', "bic")] = el(r(S), 1, 6)
		
		quietly melogit y `m`m'' || person:, asis from(B_newpersons) iterate(0)
		matrix S`m'[`i', colnumb(S`m', "newpersons")] = -2 * e(ll)
		
		quietly melogit y `m`m'' || person:, asis from(B_newitems) iterate(0)
		matrix S`m'[`i', colnumb(S`m', "newitems")] = -2 * e(ll)
		
		quietly melogit y `m`m'' || person:, asis from(B_newboth) iterate(0)
		matrix S`m'[`i', colnumb(S`m', "newboth")] = -2 * e(ll)
	
		matrix S`m'[`i', colnumb(S`m', "failures")] = `failures'
	
	}
	
}

timer off 1
timer list

forvalues m = 1/3 {
	clear
	svmat double S`m', names(col)
	generate model = `m'
	export delimited using "sim1_result_m`m'_`s'.csv", replace
}

