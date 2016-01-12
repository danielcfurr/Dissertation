// Conduct simulation ----------------------------------------------------------

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
generate newpersons = .
generate newitems = .
generate newboth = .

mkmat *, matrix(S1)
matrix define S2 = S1
matrix define S3 = S1

local m1 "w1 w2 x2-x4"
local m2 "w1 w2 x2-x4 c.x2#c.x3"
local m3 "w1 w2 x2-x4 c.x2#c.x3 c.x3#c.x4"
local N = _N

clear

forvalues i = 1/`N' {

	display "$S_TIME: Started `i'"

	set seed `=el(S1, `i', 1)'
	sim_dataset, tau(`=el(S1, `i', 2)')

	forvalues m = 1/3 {
	
		quietly melogit y_newpersons `m`m'' || person:
		matrix B_newpersons = e(b)
		local cols : colfullnames B_newpersons 
		matrix colnames B_newpersons = `=subinstr("`cols'", "_newpersons", "", .)'

		quietly melogit y_newitems `m`m'' || person:
		matrix B_newitems = e(b)
		local cols : colfullnames B_newitems 
		matrix colnames B_newitems = `=subinstr("`cols'", "_newitems", "", .)'
		
		quietly melogit y_newboth `m`m'' || person:
		matrix B_newboth = e(b)
		local cols : colfullnames B_newboth 
		matrix colnames B_newboth = `=subinstr("`cols'", "_newboth", "", .)'

		quietly melogit y `m`m'' || person:
		matrix S`m'[`i', 3] = -2 * e(ll)
		
		quietly estat ic
		matrix S`m'[`i', 4] = el(r(S), 1, 5)
		matrix S`m'[`i', 5] = el(r(S), 1, 6)
		
		quietly melogit y `m`m'' || person:, asis from(B_newpersons) iterate(0)
		matrix S`m'[`i', 6] = -2 * e(ll)
		
		quietly melogit y `m`m'' || person:, asis from(B_newitems) iterate(0)
		matrix S`m'[`i', 7] = -2 * e(ll)
		
		quietly melogit y `m`m'' || person:, asis from(B_newboth) iterate(0)
		matrix S`m'[`i', 8] = -2 * e(ll)
	
	}
	
}

forvalues m = 1/3 {
	clear
	svmat double S`m', names(col)
	generate model = `m'
	save  "sim1_`s'_`m'.dta", replace
}

