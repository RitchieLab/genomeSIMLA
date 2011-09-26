//
// C++ Interface: ese
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#define HOMOZYGOTE1 	1
#define HETEROZYGOTE 	2
#define HOMOZYGOTE2 	3
#define NODATA      	0


/*! \mainpage 
 * \section app_intro MDR-PDT Introduction
 * <P>This Version of MDR-PDT represents the first stage of a modular analysis tool for SNP data. In the current version, 
 * the tool is able to perform analysis on pedigree data based on evaluation of the T-Statistic Martin, et al (2006). 
 * 
 * \section mdrpdt_userguid User's Guide
 * <P><H2>Datasets Appropriate for MDR-PDT</H2>
 * While MDR (Ritchie, et al., 2001) was designed to discover locus-disease associations in case-control data, MDR-PDT was designed to discover
 * locus-disease associations in pedigree data using the genotype-PDT statistic (Martin et al., 2003). This version of MDR-PDT is capable of analyzing
 * input in standard ped format, which allows for a mix of trio data and discordant sibships. In both types of family structures, the informative
 * element is the transfer of genotypes to an affected child as compared with the genotypes transferred to an unaffected child. In the trio data, a
 * virtual unaffected sibling is created who is the genotypic complement of the affected child as derived from the original affected child and the two
 * parents. 
 * <P>Each locus must contain only two alleles or three genotypes and may also include a 0 indicating that the data is missing.
 * <P><H2>Pedigree Format</H2>
 * <P><H3>Trios</H3>
 * Trio families are represented by two parents and a single affected child. During the preparation of a triad, a virtual unaffected sibling is
 * produced. The production of virtual siblings is documented in the pedigree report. These entries include the description of the parents, the
 * affected sibling and the newly created virtual sibling. The new children will be given new IDs. It should be noted that missing data from any
 * parent will be propagated to both siblings.
 * <P><H3>Discordant Sibships</H3>
 * A discordant sibship is a sibship in which there one ore more affected individuals and at least one unaffected. 
 * 
 * ESE currently accepts a single command line argument: <I>the configuration file.</I>
 * <P>This file contains all the details required to execute a single run. Options include:
 * <UL>
 * <LI><B>LOCI_COUNT [integer value]</B><P> Value represents max loci count for models being evaluated.
 * <LI><B>ANALYSISSTYLE [selection]</B> 
 *     <UL><LI>NOANALYSIS - Performs no analysis</LI>
 *         <LI>MAXIMIZEDIFF - Performs analysis using the maximized difference forumala</LI></UL>
 * <LI><P><B>BESTDIFFERENCETHRESHOLD [integer value]</B><P> The value sets the threshold for which a model is considered worth reporting
 * <LI><B>VERBOSE [on/off]</B> <P>Turns the verbose mode on/off. 
 * <LI><B>INPUTFORMAT [selection]</B> 
 *     <UL><LI>MDRFORMAT - Format used in most versions of MDR (no headers at this time)</LI>
 *         <LI>INPUTSPACEDELIMITED - Format in which snps are aligned in rows (individuals are columns) separated by space</LI>
 *         <LI>INPUTLABELEDSPACEDELIMITED - Same as above except reading in a label</LI>
 *         <LI>INPUTBINARYINPUT - Binary input</LI></UL>
 * <LI><P><B>INPUTFILE [file name]</B> <P>The name of the file where the snp information can be found
 * <LI><B>REPORTLOCUSUSE [selection]</B> 
 *     <UL><LI>ECHOLOCUS - Reports the inclusion status at each snp</LI>
 *         <LI>ECHOINCLUDEDLOCUS* - Reports only loci that are included</LI>
 *         <LI>ECHOEXCLUDEDLOCUS* - Reports only the loci that are not included</LI></UL>
 * <LI><P><B>LOCUSREPORTFILE [file name]</B> <P>Specifies the name of the file where the <I>Locus Use Report</I> is stored. If no file is specified, loci report is written to std out.
 * <LI><B>REPORTMODELCOUNT [integer value]</B> <P>Value indicates how many models are to be reported. If 1 or more is specified, results
 * will be returned sorted by highest value (so it will be top <I>n</I> models listed. If 0 is specified, all models will be listed that 
 * meet or exceed the threshold. 
 * <LI><B>PERFORMSORT [yes/no]</B> <P>Indicates that the report should be sorted (only used if reportmodelcount = 0).
 * <LI><B>REPORTMODELSOUTPUTFILE [file name]</B> <P>Specify name of file to be used for reporting noteworthy models. If no file is specified, results are written to std out.
 * <P><P>*items with an asterisk beside them are not currently implimented
 * </UL>
 * \section ese_algorithm Algorithm
 * <P>
 * <UL><LI>Load into memory, all snps for each individual that are to be considered
 *     <LI>Build models that represent 1 or 2 SNPs and their difference value (/aff - unaff/) across the 
 * 		   different genotypes
 *     <LI>Perform some basic analysis on the models, attempting to reduce the number of models to 
 *         only the most interesting set
 *     <UL><LI><I>this is where things are a bit fuzzy</I> Reduce models where a single SNP occurs over and over to 
 *         just the most important pair.
 *         <LI>Remove any model that doesn't meet the configuration's threshold.
 *         <LI>Try and identify models where the pair of SNPs generally are unaffected unless they are explicitly paired
 *         together.
 *         </UL>
 *     </UL>
 * \section ese_design Design
 * A small number of classes will be used to represent the genotypic information found in the data set.
 * <P>Classes:
 * <UL><LI><B>Bitset </B>Component from the standard library which encapsulates bit manipulation tasks on
 *     predefined sized sets of bits. A single bitset will be used for each genotype (possibly two per genetotype
 *     1 affected, 1 unaffected)
 *     <LI><B>SnpAligned </B>In memory representation of the snp/model. A single snp will have N+1 bitsets to 
 *     represent each of the genotypes available (where N is the number of genotypes). Bitset 0 is keep up with
 *     individuals that were missing data for that snp/model.
 *     <LI><B>SnpRepository </B> Home for each of the SNPs. This will probably be a Blitz:Array 
 *     <LI><B>ModelRepository </B>This will contain all models and will be used to perform the various 
 *     decisions about whether to retain a model, attempt to evaluate their importance and even record which
 *     models have been discarded and why. 
 * </UL>
 * <P>Interaction:
 * <UL><LI>When a dataset is loaded, each SNP will be represented in memory as 2(N+1) bitsets representing each of the 
 *     various genotypes for different users. 
 *     <LI>To construct a model, 1 or 2 SNPs are evaluated for maximum difference. 
 *     <UL><LI>For single SNPs, this is done by simply counting the ones in the bitset.
 *         <LI>For 2 SNP models, each of the SNPs will have the common genotypes &d and the resulting bitset's 1s counted
 *     </UL>
 *     <LI>Once the model has been constructed, it will be submitted to the repository which can make decisions
 *     about the model's viability in comparison with the requirements of the user's configuration. If ESE becomes
 *     designed to run in parallel, a specialized repository could transmit viable models back to the Master
 *     <LI>Reduction can be performed as a function of the complete model repository
 *     <LI>In order to retain the current behavior of ESE, intermediate representations of the repository will be 
 *     dumped in binary format to disk upon completion of processing. This binary data can then be used later
 *     for reduction or analysis when running partial tasks on separate machines. 
 * </UL>
 */
