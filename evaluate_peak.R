'Evaluate a given m/z value string for possible molecular formulas

Usage:
  evaluate_peak.R --smiles <smiles_string> --mz_val <mz> --tolerance <tol>
  evaluate_peak.R --formula <formula_string> --mz_val <mz> --tolerance <tol>

  evaluate_peak.R (-h | --help)
  evaluate_peak.R --version

Options:
  -h --help     Show this screen.
' -> doc

###################################################################################################
suppressMessages(library(rcdk))
suppressMessages(library(docopt))
suppressMessages(library(stringr))
args <- docopt(doc)
args$tol <- as.numeric(args$tol)
args$mz <- as.numeric(args$mz)
###################################################################################################
# Initialize molecule from formula/smiles
sp <- get.smiles.parser()

if (args$smiles){
    molecule <- parse.smiles(args$smiles_string)[[1]]
    convert.implicit.to.explicit(molecule)
    formula <- get.mol2formula(molecule, charge=0)
    print(formula)
} else {
    formula <- get.formula(args$formula_string, charge=0)
}

iso <- formula@isotopes
iso_rows <- dim(iso)[1]
element_list <- list()
for (i in 1:iso_rows) {
    row <- iso[i, ]
    element <- row[1][[1]]
    count <- row[2][[1]]
    element_list[[i]] <- c(element,0,count)
}

###################################################################################################
# Generate possible formulas and print

# validation was giving me issues...
possible <- generate.formula(args$mz,window=args$tol,elements=element_list, validation = FALSE)
for (i in possible) {
    print(i)
}
