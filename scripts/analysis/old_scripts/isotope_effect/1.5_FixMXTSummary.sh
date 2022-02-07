#!/bin/bash

# intention: change the MXT2Summary file to get rid of the entries containing asterisks (*); necessary for the next script (maybe change the CreateMXTSummary script!)

cp MXT2Summary.txt MXT2Summary_backup.txt
egrep -v '[*]' MXT2Summary.txt > MXT2Summary.tmp
mv MXT2Summary.tmp MXT2Summary.txt
