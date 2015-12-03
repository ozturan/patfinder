
#Patfinder

This project consists of a simple python program that enables one to look up for designed mutations in a protein sequence. Best way to express the function of the program is through an example:

Suppose you have two proteins and you know that they are interacting. Also, you know the predicted (or intrinsic) interaction sites of these two proteins. You want investigate their relations through creating mutations in their active sites. Of course first you want to do a BLAST search to see any similar se- quence that corresponds a residue that is changed in exactly the same active site. You download several number of similar sequences and do a multiple sequence alignment (most likely CLUSTALO). And then you look for aligned sequences and search for a pattern. For example your pattern can be like this:

> GTPSTKLYGDVNDDGKVNSTDAVALKRYVLRSGISINTDNADLNE
> DGRVNSTDLGILKRYILKEIDTLPYKN

You see that in the query sequence **ST-KR** pattern is existed. These pattern corresponds to an active site and If you want to do research on this site in terms binding activity, you need to make mutations. In other words you need to find meaningful mutations based on similar sequences. And that’s where PatFinder comes in.
