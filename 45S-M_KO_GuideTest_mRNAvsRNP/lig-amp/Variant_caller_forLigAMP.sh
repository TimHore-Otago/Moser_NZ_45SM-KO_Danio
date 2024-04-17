#!/bin/bash

# execute by first giving permission and then running it as a nohup: 
# chmod +x Variant_caller_forLigAMP.sh &&
# nohup ./Variant_caller_forLigAMP.sh



# # # # # # # # # # # # # # # # # # # # # #  # # # # # # # # # # # # # # # # # # # # #
#                                                                                    #
#                   Shell Version which requires to be located in the WDir           #
#                                                                                    #
####                                   part00                                     ####
###                                                                                ###
##         !         Make sure Variant_lister_V8.sh is in working dir       !       ##
#                                                                                    #
# # # # # # # # # # #  just to be save do it in several steps.   # # # # # # # # # # # 

## Theory:
    #   This version only works with lig amp reads where the indexing primers are randomly assigned to either end of the amplicon.
    #   The pipeline is essentially the same as for the lange handle indexed samples except for initial steps which sort R1/R2 into correct files based on primer not ilumina label.
    #   Cosequently the R1/R2 (i.e. 1:N:0/2:N:0) tag are applied to the specific read name. 
    #   
    #
    #
    #   on-target.  (quite lenient allowing for 3 mismatches and any 31kmer match to) 
    #       GGCAGGTCAGCGATAACCAGTGACGGCTGCCCAGACCTCCGGGAGAGAGAGACCCAAGAGGTCACGTCCCGCCCGGCTGTG
    #       AACCTGAGCTCGGGGGCCATGGATAGTCAGGGTTAGCGGCCGCTTCTTGCTTCTCCCTTGCCGAATCGGGTGCCGATGGTG
    #       ACCCGCGACTGACTTTGCCTTCGGCAGCCGGAGCGGTGCTGCTGGTCACCCACAGCTGACCTTTCTGGCTGGCACGCGGGG
    #       TGGCACTTCTCCCAAGCCTCAGGCAACCAGGCAAGCCAAAGCTCCACAAGCCCAA
    # 
    # Indel Variants: 
    #   I.  If sgRNA#1 and sgRNA#4 loci are indel free -> perfect'ish
    #       GTCACGTCCCGCCCGGCTGT,TTCTGGCTGGCACGCGGGGT
    #       (overruled by presence on other lists)
    #   II. If R1 mod + contains pos90-151.     R2 NA.  
    #       GAACCTGAGCTCGGGGGCCATGGATAGTCAGGGTTAGCGGCCGCTTCTTGCTTCTCCCTTG
    #   V.  If R1 contains  >pos167-229.         R2 NA.
    #       ATGGTGACCCGCGACTGACTTTGCCTTCGGCAGCCGGAGCGGTGCTGCTGGTCACCCACAGCT
    #   III.If perfect but R 1 contains R2 tail and vice versa. Not always detectable because of the 17bp gap!
    #       
    #   IV. If R2 mod + contains pos166-229.    R1 NA.
    #       ATGGTGACCCGCGACTGACTTTGCCTTCGGCAGCCGGAGCGGTGCTGCTGGTCACCCACAGCT
    #   VI. If R2  pos84-151.                   R1 NA.
    #       CGGCTGTGAACCTGAGCTCGGGGGCCATGGATAGTCAGGGTTAGCGGCCGCTTCTTGCTTCTCCCTTG        
    #   VII.If R1 contains  pos245-296
    #       ACGCGGGGTGGCACTTCTCCCAAGCCTCAGGCAACCAGGCAAGCCAAAGCT
    #       If R2 contains  pos21-80
    #       CGATAACCAGTGACGGCTGCCCAGACCTCCGGGAGAGAGAGACCCAAGAGGTCACGTCC

# NFO # # NFO # # NFO # # NFO # # NFO # # NFO # 
#casVar_caller_V6 for Lig-amp reads
echo "- - - - - - - " >/dev/tty

echo "! ! ! Warning ! ! ! " >/dev/tty
echo "- - - - - - - " >/dev/tty
echo "Running this script might take up to several hours.
Substantial storage capacity is required.
To reduce system load, the number of parallel processes triggered by this script are restricted. " >/dev/tty
echo "- - - - - - - " >/dev/tty
echo "! ! ! Warning ! ! ! " >/dev/tty
echo "- - - - - - - " >/dev/tty

while true; do
    # Prompt the user to confirm script execution
    echo " Do you wish to continue? (Yes/No)" >/dev/tty
    read execute_choice >/dev/tty

    case "$execute_choice" in
        [yY]|[yY][eE][sS])
            # If the user chooses Yes or yes, break out of the loop and execute the script
            break
            ;;
        [nN]|[nN][oO])
            # If the user chooses No or no, exit the script
            echo "Shell script execution canceled." >/dev/tty
            exit 0
            ;;
        *)
            # If the input is invalid, prompt the user again
            echo "Invalid choice. Please type 'Yes' or 'No'." >/dev/tty
            ;;
    esac
done

date &&
echo "- - - - -"  > /dev/tty &&
echo "Commencing Cas9 Indel Variant Identification (8 Steps) " > /dev/tty &&
echo "- - - - -"  > /dev/tty &&
date > /dev/tty &&
echo "- - - - -"  > /dev/tty &&
echo "Step 1/8 - unzipping" > /dev/tty &&
mkdir ./temp &&
mkdir ./output &&
mv ./raw_seq_files/*gz ./temp/ &&

cd ./temp/ &&

gunzip *.gz &&

#if [[ -n $(ls *_001.fastq 2>/dev/null) ]]; then
#    for file in *_001.fastq; do
#        mv -- "$file" "${file%_001.fastq}.fastq"
#    done
#fi
wait
echo "Step 2/8 - reassigning R1/R2 status " > /dev/tty &&

# concatenate
for file1 in *R1_001.fastq; do
    file2="${file1/R1_001/R2_001}"
    
    diff_count=$(diff <(echo "$file1") <(echo "$file2") | grep -c '[><]')

    if [ "$diff_count" -eq 2 ]; then
        merged_filename=$(echo "$file1" | sed -E "s/^(.*)(.)_R1_001/\1\2_conc/")
        cat "$file1" "$file2" > "$merged_filename"
        rm "$file1" "$file2"
    fi
done
wait
# re-assign R1/R2 (requires the rcomp flag to be false or else max deletion reads will be misslabelled)
#R1
max_processes=6
running_processes=0
for filename in ./*_conc.fastq; do
    name=$(basename ${filename} .fastq)
    /Applications/bbmap/bbduk.sh in=$filename outm=./${name}_onR1.fastq literal=AAACCTCGACGGCAGGTCAG k=20 hdist=2 rieb=f rcomp=f&
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done
wait 
#R2
max_processes=6
running_processes=0
for filename in ./*_conc.fastq; do
    name=$(basename ${filename} .fastq)
    /Applications/bbmap/bbduk.sh in=$filename outm=./${name}_onR2.fastq literal=GCTTAGTTGTTGGGCTTGTGG k=21 hdist=2 rieb=f rcomp=f&
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done
wait 

## re-labelling read names correctly
for file in *_onR1.fastq; do
    sed 's/2:N:0/1:N:0/g' "$file" > "$file.tmp" && mv "$file.tmp" "$file"
done
wait 
for file in *_onR2.fastq; do
    sed 's/1:N:0/2:N:0/g' "$file" > "$file.tmp" && mv "$file.tmp" "$file"
done
wait 

# re-pairing using repair.sh ( since the process above may produce singletons)
# somehow as a consequence of the lig amp procedure there are a lot of read pairs where one is perfect and the other is a mess. thus, lots of singles. 
# while it sucks to loose those reads and the infromation they carry this needs to happen for the remaining pipeline to work!
max_processes=6
running_processes=0
for filename1 in ./*_onR1.fastq; do
    name=$(basename ${filename1} _onR1.fastq)
    filename2="${name}_onR2.fastq"
    if [[ -e $filename2 ]]; then
        /Applications/bbmap/repair.sh in1=$filename1 in2=$filename2 out=./${name}_rep.fastq   repair=t& #outs=./${name}_Singles.fastq -> add if you want to see the rejected singles  
        ((running_processes++))
    fi
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done 
wait


#Trimming ( pairing was done in the above repair step)
max_processes=6
running_processes=0
for filename1 in ./*_rep.fastq; do
    name=$(basename ${filename1} _rep.fastq)
    /Applications/bbmap/bbduk.sh in=$filename1 out=./${name}_pr.fq ref=/Applications/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 minlen=111 &
    ((running_processes++))
    
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done 
wait
echo "Step 3/8 - filtering for on-target read pairs" > /dev/tty &&
#on-target filtering  -> *_on.fq  # somewhat redundant in this pipeline since the primer filters in the read assignment will have covered this(mostly)
max_processes=6
running_processes=0
for filename in ./*pr.fq; do
    name=$(basename ${filename} .fq)
    /Applications/bbmap/bbduk.sh in=$filename outm=./${name}_on.fq literal=GGCAGGTCAGCGATAACCAGTGACGGCTGCCCAGACCTCCGGGAGAGAGAGACCCAAGAGGTCACGTCCCGCCCGGCTGTGAACCTGAGCTCGGGGGCCATGGATAGTCAGGGTTAGCGGCCGCTTCTTGCTTCTCCCTTGCCGAATCGGGTGCCGATGGTGACCCGCGACTGACTTTGCCTTCGGCAGCCGGAGCGGTGCTGCTGGTCACCCACAGCTGACCTTTCTGGCTGGCACGCGGGGTGGCACTTCTCCCAAGCCTCAGGCAACCAGGCAAGCCAAAGCTCCACAAGCCCAA k=31 hdist=3 qin=33 rieb=f &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done 
wait
echo "Step 4/8 - seperating on-target R1/R2 reads" > /dev/tty &&
###split all onT R1/R2
max_processes=6
running_processes=0
for filename1 in ./*_on.fq; do
    name=$(basename ${filename1} _on.fq)
    /Applications/bbmap/bbduk.sh in=$filename1 out1=./${name}_onR1.fq out2=./${name}_onR2.fq &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done 
wait
echo "Step 5/8 - assigning all indel variants possibly occuring in the respective read" > /dev/tty &&
### find I (i.e perfectish)
max_processes=6
running_processes=0
for filename in ./*_on.fq; do
    name=$(basename ${filename} .fq)
    /Applications/bbmap/bbduk.sh in=$filename outm=./${name}_perfish.fq literal=GTCACGTCCCGCCCGGCTGT,TTCTGGCTGGCACGCGGGGT k=20 hdist=0 qin=33 rieb=f &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done 
wait




### find II  (R1only)
#find mod R1
max_processes=6
running_processes=0
for filename in ./*_onR1.fq; do
    name=$(basename ${filename} _onR1.fq)
    /Applications/bbmap/bbduk.sh in=$filename out=./${name}_onR1mod.fq literal=GTCACGTCCCGCCCGGCTGT k=20 hdist=0 &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done
wait 
#contains pos99-151.
max_processes=6
running_processes=0
for filename1 in ./*_onR1mod.fq; do
    name=$(basename ${filename1} _onR1mod.fq)
    /Applications/bbmap/bbduk.sh in1=$filename1 outm=./${name}_sg1indel.fq literal=GAACCTGAGCTCGGGGGCCATGGATAGTCAGGGTTAGCGGCCGCTTCTTGCTTCTCCCTTG k=23 hdist=1 &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done
wait 
### find IV (R2only)
#find mod R2
max_processes=6
running_processes=0
for filename in ./*_onR2.fq; do
    name=$(basename ${filename} _onR2.fq)
    /Applications/bbmap/bbduk.sh in=$filename out=./${name}_onR2mod.fq literal=TTCTGGCTGGCACGCGGGGT k=20 hdist=0 qin=33 &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done
wait 
#contains pos190-229. 
max_processes=6
running_processes=0
for filename1 in ./*_onR2mod.fq; do
    name=$(basename ${filename1} _onR2mod.fq)
    /Applications/bbmap/bbduk.sh in1=$filename1 outm=./${name}_sg4indel.fq literal=ATGGTGACCCGCGACTGACTTTGCCTTCGGCAGCCGGAGCGGTGCTGCTGGTCACCCACAGCT k=23 hdist=1 &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done
wait 
###find V  (R1only) 
#too long-> will not flag short reads with V
#too short-> will possibly flag II, III as V
max_processes=6
running_processes=0
for filename1 in ./*_onR1.fq; do
    name=$(basename ${filename1} _onR1.fq)
    /Applications/bbmap/bbduk.sh in1=$filename1 outm=./${name}_ldel1to3.fq literal=ATGGTGACCCGCGACTGACTTTGCCTTCGGCAGCCGGAGCGGTGCTGCTGGTCACCCACAGCT k=25 hdist=2 &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done
wait 
###find VI  (R2only)
max_processes=6
running_processes=0
for filename1 in ./*_onR2.fq; do
    name=$(basename ${filename1} _onR2.fq)
    /Applications/bbmap/bbduk.sh in1=$filename1 outm=./${name}_ldel3to4.fq literal=CGGCTGTGAACCTGAGCTCGGGGGCCATGGATAGTCAGGGTTAGCGGCCGCTTCTTGCTTCTCCCTTG k=25 hdist=2 &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done
wait 
###find VII 
#VIIa  (R1)
max_processes=6
running_processes=0
for filename1 in ./*_onR1.fq; do
    name=$(basename ${filename1} _onR1.fq)
    /Applications/bbmap/bbduk.sh in1=$filename1 outm=./${name}_ldel1to4R1.fq literal=ACGCGGGGTGGCACTTCTCCCAAGCCTCAGGCAACCAGGCAAGCCAAAGCT k=25 hdist=2 &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done
wait 
#VIIb  (R2)
max_processes=6
running_processes=0
for filename1 in ./*_onR2.fq; do
    name=$(basename ${filename1} _onR2.fq)
    /Applications/bbmap/bbduk.sh in1=$filename1 outm=./${name}_ldel1to4R2.fq literal=CGATAACCAGTGACGGCTGCCCAGACCTCCGGGAGAGAGAGACCCAAGAGGTCACGTCC k=25 hdist=2 &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done
wait
### find III (R2 no ldel3-4)
max_processes=6
running_processes=0
for filename in ./*_onR2.fq; do
    name=$(basename ${filename} _onR2.fq)
    /Applications/bbmap/bbduk.sh in=$filename outm=./${name}_R2noVI.fq literal=GCGACTGACTTTGCCTTCGGCAGCCGGAGCGGTGCTGCTG k=25 hdist=2 &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done
wait 
max_processes=6
running_processes=0
for filename in ./*_R2noVI.fq; do
    name=$(basename ${filename} _R2noVI.fq)
    /Applications/bbmap/bbduk.sh in=$filename outm=./${name}_R2wIII.fq literal=CCGCTTCTTGCTTCTCCCTTGCCGAATCGG k=11 hdist=1 &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done
wait 
### find III R1 no ldel3-4 
max_processes=6
running_processes=0
for filename in ./*_onR1.fq; do
    name=$(basename ${filename} _onR1.fq)
    /Applications/bbmap/bbduk.sh in=$filename outm=./${name}_R1noV.fq literal=GGATAGTCAGGGTTAGCGGCCGCTTCTTGCTTCTCCCTTG k=25 hdist=2 &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done
wait 
max_processes=6
running_processes=0
for filename in ./*_R1noV.fq; do
    name=$(basename ${filename} _R1noV.fq)
    /Applications/bbmap/bbduk.sh in=$filename outm=./${name}_R1wIII.fq literal=CCGATGGTGACCCGCGACTGACTTTGCCTT  k=11 hdist=1 &
    ((running_processes++))
    if [[ $running_processes -ge $max_processes ]]; then
        wait
        running_processes=0
    fi
done
wait

#### This concludes part 1 ####
echo "Step 6/8 - making lists with reads assigned to each indel variants" > /dev/tty &&
# # # # # # # # # # # # # # # # # # # # # #  # # # # # # # # # # # # # # # # # # # # #
#                                                                                    #  
####                                   part02                                     ####
#                                                                                    #
# # # # # # # # # # #  just to be save do it in several steps.   # # # # # # # # # # # 


#get readlists
for file in *.fq    
do
    grep "^@" "$file" > "${file%.fq}_rl.txt" &    
done 
wait
for file in *_ldel1to3_rl.txt; 
do mv "$file" "${file/_ldel1to3_rl/_ldel1to3_V_rl}" 
done 
wait
for file in *_ldel1to4R1_rl.txt; 
do mv "$file" "${file/_ldel1to4R1_rl/_ldel1to4R1_VIIr1_rl}" 
done 
wait
for file in *_ldel1to4R2_rl.txt;
do mv "$file" "${file/_ldel1to4R2_rl/_ldel1to4R2_VIIr2_rl}" 
done 
wait
for file in *_ldel3to4_rl.txt; 
do mv "$file" "${file/_ldel3to4_rl/_ldel3to4_VI_rl}" 
done
wait
for file in *_on_perfish_rl.txt;
do mv "$file" "${file/_on_perfish_rl/_on_perfish_I_rl}" 
done 
wait
for file in *_sg1indel_rl.txt; 
do mv "$file" "${file/_sg1indel_rl/_sg1indel_II_rl}" 
done 
wait 
for file in *_sg4indel_rl.txt;
do mv "$file" "${file/_sg4indel_rl/_sg4indel_IV_rl}"
done 
wait
#just to make the filenames shorter
for file in *.txt;
do mv "$file" "${file/_S*cc_/_}"
done 
wait


#### This concludes part 2 ####

# # # # # # # # # # # # # # # # # # # # # #  # # # # # # # # # # # # # # # # # # # # #
#                                                                                    #  
####                                   part03                                     ####
#                                                                                    #
# # # # # # # # # # #  just to be save do it in several steps.   # # # # # # # # # # # 
echo "Step 7/8 - SLOW - making a presence/absence master list for all reads, based on the indel variant lists"  > /dev/tty &&
#run comparison shell (which is still unable to run as for loop) 
#give permission to execute:
cd ../ &&
wait
chmod +x Variant_lister_V8.sh &&
#execute comparison shell (content listed at bottom)
nohup ./Variant_lister_V8.sh &
wait
##  cleanup
cd ./temp &&
mkdir ./concatandpaired;mv *conc.fastq $_ &&
mv *conc_rep.fastq $_ &&
mkdir ./cas9_I;mv *fish.fq $_ &&
mkdir ./cas9_II;mv *sg1indel.fq $_ &&
mkdir ./cas9_III;mv *III.fq $_ &&
rm *_R2noVI.fq; rm *_R1noV.fq &&
mkdir ./cas9_IV;mv *sg4indel.fq $_ &&
mkdir ./cas9_V;mv *ldel1to3.fq $_ &&
mkdir ./cas9_VI;mv *ldel3to4.fq $_ &&
mkdir ./cas9_VIIa;mv *ldel1to4R1.fq $_ &&
mkdir ./cas9_VIIb;mv *ldel1to4R2.fq $_ &&
mkdir ./R1on;mv *onR1.fq $_ ;mv *onR1mod.fq $_ &&
mkdir ./R2on;mv *onR2.fq $_;mv *onR2mod.fq $_ &&
mkdir ./onTarget;mv *on.fq $_ &&
mkdir ./lists;mv *rl.txt $_ &&
mkdir ./TrimPaired;mv *pr.fq $_
mv *R1.fastq ../raw_seq_files ;mv *R2.fastq ../raw_seq_files &&
mv *put.txt ../output &&
wait &&

cd ../

wait &&
gzip ./temp/*/*.fq &&
wait &&
gzip ./temp/*/*.txt &&
wait &&
gzip ./raw_seq_files/*.fastq &&
wait &&
date &&

echo "- - - - -"  > /dev/tty &&
date > /dev/tty &&
echo "- - - - -"  > /dev/tty 
echo "All Done!" > /dev/tty &&
echo "- - - - -"  > /dev/tty 
echo "The 'temp' folder is no longer required and can be deleted." > /dev/tty &&
#!/bin/bash

while true; do
    # Prompt the user
    echo "Do you want to delete the 'temp' folder? (Yes/No)" >/dev/tty
    read choice >/dev/tty

    # Check the user's choice
    case "$choice" in
        [yY]|[yY][eE][sS])
            # If the user types Yes or yes, delete the folder
            rm -r ./temp
            echo "Folder 'temp' has been deleted.">/dev/tty
            break  # Exit the loop
            ;;
        [nN]|[nN][oO])
            # If the user types No or no, do nothing
            echo "Folder 'temp' will not be deleted.">/dev/tty
            break  # Exit the loop
            ;;
        *)
            # If the input doesn't match Yes or No, display an error message
            echo "Invalid choice. Please type 'Yes' or 'No'." >/dev/tty
            ;;
    esac
done



##################################### all done ###################################