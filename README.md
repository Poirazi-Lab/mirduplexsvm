# README

(migrated repo from Karathanasis N, personal repo: [bitbucket.org/nestoras_karathanasis/mirduplexsvm/](https://bitbucket.org/nestoras_karathanasis/mirduplexsvm/))

## What is this repository for?

Provide MiRduplexSVM code to the scientific community.

Karathanasis N, Tsamardinos I, Poirazi P (2015) MiRduplexSVM: A High-Performing MiRNA-Duplex Prediction and Evaluation Methodology. PLOS ONE 10(5): e0126151. doi: [10.1371/journal.pone.0126151](https://doi.org/10.1371/journal.pone.0126151)

## Paper Summary
We address the problem of predicting the position of a miRNA duplex on a microRNA hairpin via the development and application of a novel SVM-based methodology. Our method combines a unique problem representation and an unbiased optimization protocol to learn from mirBase19.0 an accurate predictive model, termed MiRduplexSVM. This is the first model that provides precise information about all four ends of the miRNA duplex. We show that (a) our method outperforms four state-of-the-art tools, namely MaturePred, MiRPara, MatureBayes, MiRdup as well as a Simple Geometric Locator when applied on the same training datasets employed for each tool and evaluated on a common blind test set. (b) In all comparisons, MiRduplexSVM shows superior performance, achieving up to a 60% increase in prediction accuracy for mammalian hairpins and can generalize very well on plant hairpins, without any special optimization. (c) The tool has a number of important applications such as the ability to accurately predict the miRNA or the miRNA*, given the opposite strand of a duplex. Its performance on this task is superior to the 2nts overhang rule commonly used in computational studies and similar to that of a comparative genomic approach, without the need for prior knowledge or the complexity of performing multiple alignments. Finally, it is able to evaluate novel, potential miRNAs found either computationally or experimentally. In relation with recent confidence evaluation methods used in miRBase, MiRduplexSVM was successful in identifying high confidence potential miRNAs.

* Version 1

## How do I get set up?

- Download the code.
- No Configuration is needed
- Dependencies: RNAfold program from Vienna RNA Package, version 1.8.5. LibSVM, version 3.1
- Follow the [instructions](https://github.com/Poirazi-Lab/mirduplexsvm/blob/main/instructions.md)

## Who do I talk to?

- Admin
