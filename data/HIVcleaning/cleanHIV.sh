#!/bin/bash
# initial cleaning of paediatric HIV/ART data

# HIV file from http://aidsinfo.unaids.org, only numbers available for prevalence of child infetions
sed '1d' KHIVraw.txt > tmp/K2.txt  # ditch top line
awk -F", *" '!/[\.\.\.]/' tmp/K2.txt > tmp/K3a.txt #ditch anything with the ... field
awk -F", *" '!/Kyrgyzstan/' tmp/K3a.txt > tmp/K3.txt #ditch Kyrgyzstan which has uninterpretable data
cut -d ',' -f1 tmp/K3.txt > tmp/K4a.txt 	#print first column to file
cut -d ',' -f2- tmp/K3.txt | tr -d ' ' > tmp/K4b.txt #ditch all spaces beyond 1st col
cat tmp/K4b.txt | sed s/'PMTCTnumeratorfor2015wasunavailableatthetimeofdevelopmentofprojection,*'//g > tmp/K4b1.txt #ditch these extra fields
cat tmp/K4b1.txt | sed s/'Antiretroviraltherapydatawasnotavailableatthetimeofpublication,*'//g > tmp/K4b2.txt #ditch these extra fields
paste -d ',' tmp/K4a.txt tmp/K4b2.txt > tmp/K5.txt # join the two previous files
cp tmp/K5.txt KHcleaner.csv


# record missing countries
awk -F", *" '/[\.\.\.]/' tmp/K2.txt > tmp/K3arem.txt #lines ditched
cut -d ',' -f1 tmp/K3arem.txt > tmp/countrydrops.txt 	#print first column to file

# ART coverage file from http://aidsinfo.unaids.org
sed '1d' KARTraw.csv > tmp/KA2.txt  # ditch top line
awk -F", *" '!/[\.\.\.]/' tmp/KA2.txt > tmp/KA3a.txt #ditch anything with the ... field
awk -F", *" '!/Kyrgyzstan/' tmp/KA3a.txt > tmp/KA3.txt #ditch Kyrgyzstan which has uninterpretable data
cut -d ',' -f1 tmp/KA3.txt > tmp/KA4a.txt 	#print first column to file
cut -d ',' -f2- tmp/KA3.txt | tr -d ' ' > tmp/KA4b.txt #ditch all spaces beyond 1st col
cat tmp/KA4b.txt | sed s/'PMTCTnumeratorfor2015wasunavailableatthetimeofdevelopmentofprojection,*'//g > tmp/KA4b1.txt #ditch these extra fields
paste -d ',' tmp/KA4a.txt tmp/KA4b1.txt > tmp/KA5.txt # join the two previous files
cp tmp/KA5.txt KAcleaner.csv

