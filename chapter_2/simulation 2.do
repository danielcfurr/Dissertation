local s = 1
run sim_dataset.do

import delimited "random.csv", varnames(nonames) clear
rename v`s' seed
drop v*
generate double tau = 0
replace tau = .1 in 2
replace tau = .5 in 3
replace tau = 1 in 4
fillin tau seed
drop _fillin
generate insample = .
generate aic = .
generate bic = .
generate newitems = .

mkmat *, matrix(S1)
matrix define S2 = S1
matrix define S3 = S1

local m1 "x2-x4"
local m2 "x2-x4 c.x2#c.x3"
local m3 "x2-x4 c.x2#c.x3 c.x3#c.x4"
local N = _N

timer clear
timer on 1

forvalues i = 1/`N' {

	display "$S_TIME: Started `i'"

	set seed `=el(S1, `i', 1)'
	local retry = 1
	local ntries = 0

	while `retry' == 1 {

		sim_dataset, tau(`=el(S1, `i', 2)')
		
		local retry = 0
		local ++ntries
		quietly tabulate item
		local I = r(r)
		
		preserve

			collapse (sum) y y_newitems, by(person)

			// Check that sum scores are not all correct or all incorrect
			quietly summarize y
			if (r(min) == 0 | r(max) == `I') local retry = 1
			quietly summarize y_newitems
			if (r(min) == 0 | r(max) == `I') local retry = 1
			
			// Check that all other sum scores are represented
			quietly tabulate y
			if (r(r) != `I' - 1) local retry = 1
			quietly tabulate y_newitems
			if (r(r) != `I' - 1) local retry = 1
			
			//quietly levelsof y
			//local y_vals = r(levels)
			//if (`retry' == 0) {
			//	foreach lname of local y_vals {
			//		quietly count if y_newitems == `lname'
			//		if (r(N) == 0) local retry = 1
			//	}
			//}
			
			//quietly levelsof y_newitems
			//local y_vals = r(levels)
			//if (`retry' == 0) {
			//	foreach lname of local y_vals {
			//		quietly count if y == `lname'
			//		if (r(N) == 0) local retry = 1
			//	}
			//}
			
		restore

		egen sumscore = total(y), by(person)
		egen sumscore_newitems = total(y_newitems), by(person)

		forvalues m = 1/3 {

			if (`retry' == 0) {
			
				quietly melogit y_newitems i.sumscore_newitems `m`m'' || item:, iterate(100)
				if (e(converged) != 1) local retry = 1
				matrix B_newitems = e(b)
				local cols : colfullnames B_newitems 
				matrix colnames B_newitems = `=subinstr("`cols'", "_newitems", "", .)'

				quietly melogit y i.sumscore `m`m'' || item:, iterate(100)
				if (e(converged) != 1) local retry = 1
				matrix S`m'[`i', 3] = -2 * e(ll)
				
				quietly estat ic
				matrix S`m'[`i', 4] = el(r(S), 1, 5)
				matrix S`m'[`i', 5] = el(r(S), 1, 6)
				
				quietly melogit y i.sumscore `m`m'' || item:, asis from(B_newitems) iterate(0)
				matrix S`m'[`i', 6] = -2 * e(ll)

				if (`retry' == 1) display "  convergence failure--retrying"
				
			}
			
		}	
		
	}

	display "  " `ntries' " attempts at making data"
	
}

timer off 1
timer list

forvalues m = 1/3 {
	clear
	svmat double S`m', names(col)
	generate model = `m'
	save  "sim2_`s'_`m'.dta", replace
}

