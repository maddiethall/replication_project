# The relationship of coping style and social support variation to glucocorticoid metabolites in wild olive baboons (Papio anubis)

Alexander J. Pritchard, Erin R. Vogel, Rosemary A. Blersch, Ryne A. Palombit

Data were collected from two groups of in situ (wild) olive baboons, in Laikipia, Kenya. An experimental paradigm was used to collect coping style scores over a 17-month period; behavioral observations were collected with focal follows and ad libitum sampling, over a 17-month period; fecal samples were collected over a nine month period; daily rainfall and min/max temperature were collected each day that fecal sampling was ongoing. Shannon Weiner's Diversity Indices as a metric of social support calculated based on grooming dynamics across the relevant project period. Rank was calculated from social displacements. Rainfall and temperature were averaged across the two days, relative to each fecal sample's collection. Fecal sample hormones were extracted and assayed using RIAs. More details are in the associated manuscript.

## Description of the data and file structure

There are multiple columns, each line is a unique fecal sample. Some items are tied to the fecal sample (e.g., Metabolites.Conc, MAXAvg), while other values are tied to the IDs themselves (e.g., Coping, SWDIgr). Columns are: (untitled first column) & X = sequential row identifiers; id = subject identifier (initials); MAXAvg = maximum daily temperature in Celsius, averaged for the 2 days prior to fecal sample collection;	MINAvg = minimum daily temperature in Celsius, averaged for the 2 days prior to fecal sample collection; RainAvg = daily rainfall in mm, averaged for the 2 days prior to fecal sample collection; Date = Date of sample collection in YYYY-MM-DD; Time = time of sample collection in HH:MM:SS AM/PM;	RunDate = Run date of sample in DD-Mon.; Metabolites.Conc = fecal glucocorticoid metabolite concentration in a fecal sample, adjusted by weight (ng/g); RankProp = proportional social dominance rank of ID; Sex = sex of the animal subject; Group.1 = social group of the animal subject (K = Kati Kati, S = Shire); Coping = coping style score of the animal subject; SWDIgr = Shannon Weiner's diversity index for grooming of the subject; Strg = social grooming network out-strength; Deg = social grooming network out-degree.

NAs signify data that are 'not available' as they were not sampled.

There is a single *.csv file with the dataset used, in addition to R code for the relevant statistical models.

## Sharing/Access information

Please cite the original source of these data.

## Code/Software

See the R file for further details regarding analysis and procedure

Accessible using the R programming language and a *.csv reader, such as Excel or LibreOffice Calc.