/*
 * Calculate R^2 between host and symbiont variants for an S&R experiment.
 *
 * intergenomic_ld <symbiont frequencies> <symbiont VCF> <host VCF>
 *
 * Assumes all biallelic variants, and that symbiont is haploid.
 *
 * g++ -Wall -O3 --std=c++11 -o intergenomic_ld intergenomic_ld.cpp -I/usr/include/htslib -lhts -lm
 *
 */

#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <hts.h>
#include <vcf.h>

using std::cout;
using std::cerr;
using std::endl;
using std::getline;
using std::ifstream;
using std::make_pair;
using std::map;
using std::pair;
using std::string;
using std::stringstream;
using std::vector;

float calcR2(float AB, float A, float B) {
    float D = AB - (A * B);
    float r = D / sqrt(A * (1 - A) * B * (1 - B));
    return pow(r, 2);
}


int main(int argc, char * argv[]) {

    // Step 1: Read table of strain frequencies
    string strainSampleString = "";
    string hostSampleString = "";
    string host, strain;
    map<pair<string, string>, float> strainFrequencies;
    ifstream freqFile(argv[1]);
    string line, item;
    vector<string> fields;
    vector<string> fieldnames;
    int nl = -1;
    while(getline(freqFile, line)) {
        nl++;
        stringstream ss(line);
        while(getline(ss, item, '\t')) {
            fields.push_back(item);
        }
        int n = -1;
        for(const auto& f:fields) {
            n++;
            if(nl == 0) {
                if(n == 0) continue;
                fieldnames.push_back(f);
                if(n > 1) hostSampleString = hostSampleString.append(",");
                hostSampleString = hostSampleString.append(f);
            } else {
                if(n == 0) {
                    if(nl > 1) strainSampleString = strainSampleString.append(",");
                    strainSampleString = strainSampleString.append(f);
                    strain = f;
                } else {
                    host = fieldnames[n];
                    strainFrequencies[make_pair(strain, host)] = (float)atof(f.c_str());
                }
            }
        }
        fields.clear();
    }
    freqFile.close();


    /*****************************************************************/

    // Step 2: Read and store all symbiont variants
    vcfFile* symbiontVCF = bcf_open(argv[2], "r");
    bcf_hdr_t* symbiontHdr = bcf_hdr_read(symbiontVCF);
    bcf_hdr_set_samples(symbiontHdr, strainSampleString.c_str(), 0); // Only get the genotypes for particular strains
    bcf1_t* record = bcf_init();

    vector<map<string, int>> symbiontGenotypes;     // Genotype states
    vector<string> symbiontRs;                      // Variant IDs
    int i, j, ngt, nsmpl = bcf_hdr_nsamples(symbiontHdr);   // Setup for reading genotypes
    int32_t *gt_arr = NULL, ngt_arr = 0;
    int max_ploidy = 1; // Assume haploid symbiont

    // Keep vector of sample names in same order as variant file
    vector<string> symbiontSamples;
    for(int k = 0; k < nsmpl; k++) {
        symbiontSamples.push_back(string(symbiontHdr->samples[k]));
    }

    // Loop over lines in the variant file
    while(bcf_read(symbiontVCF, symbiontHdr, record) == 0) {

        // Have to tell htslib to parse the information
        bcf_unpack(record, BCF_UN_ALL);

        // Keep ordered list of variant IDs
        symbiontRs.push_back(string(record->d.id));

        // Store genotype states; key = strain, value = genotype index
        map<string, int> genotypes;

        ngt = bcf_get_genotypes(symbiontHdr, record, &gt_arr, &ngt_arr);

        for (i=0; i<nsmpl; i++) {
            int32_t *ptr = gt_arr + i*max_ploidy;
            for (j=0; j<max_ploidy; j++) {

                // if true, the sample has smaller ploidy
                if ( ptr[j]==bcf_int32_vector_end ) break;

                // missing allele: do not enter a value
                if ( bcf_gt_is_missing(ptr[j]) ) continue;

                // the VCF 0-based allele index stored in genotypes
                int allele_index = bcf_gt_allele(ptr[j]);
                if(allele_index > 1) {
                    cerr << "ERROR: Found genotype > 1 (multiallelic?) " <<
                        record->d.id << endl;
                    exit(1);
                }
                genotypes[symbiontSamples[i]] = allele_index;
            }
        }
        symbiontGenotypes.push_back(genotypes);

    }

    // Clean up
    bcf_hdr_destroy(symbiontHdr);
    bcf_close(symbiontVCF);



    /*****************************************************************/

    // Step 3: Loop through host variants and calculate LD for all combinations
    //         of host and symbiont variant

    cout << "symbiont_variant\thost_variant\tr2" << endl;

    vcfFile* hostVCF = bcf_open(argv[3], "r");
    bcf_hdr_t* hostHdr = bcf_hdr_read(hostVCF);
    bcf_hdr_set_samples(hostHdr, hostSampleString.c_str(), 0); // Only get the genotypes for particular strains

    i, j, ngt, nsmpl = bcf_hdr_nsamples(hostHdr);   // Setup for reading genotypes
    *gt_arr = NULL, ngt_arr = 0;
    max_ploidy = ngt/nsmpl;

    // Keep vector of sample names in same order as variant file
    vector<string> hostSamples;
    for(int k = 0; k < nsmpl; k++) {
        hostSamples.push_back(string(hostHdr->samples[k]));
    }

    // Loop over lines in the variant file
    map<string, vector<int>> hostGenotypes;
    while(bcf_read(hostVCF, hostHdr, record) == 0) {

        // Have to tell htslib to parse the information
        bcf_unpack(record, BCF_UN_ALL);

        // Keep ordered list of variant IDs
        string hostRs = string(record->d.id);

        // Store genotype states; key = strain, value = genotype index

        ngt = bcf_get_genotypes(hostHdr, record, &gt_arr, &ngt_arr);

        for (i=0; i<nsmpl; i++) {
            int32_t *ptr = gt_arr + i*max_ploidy;
            for (j=0; j<max_ploidy; j++) {
                // if true, the sample has smaller ploidy
                if ( ptr[j]==bcf_int32_vector_end ) break;

                // missing allele: do not enter a value
                if ( bcf_gt_is_missing(ptr[j]) ) continue;

                // the VCF 0-based allele index stored in genotypes
                int allele_index = bcf_gt_allele(ptr[j]);
                if(allele_index > 1) {
                    cerr << "ERROR: Found genotype > 1 (multiallelic?) " <<
                        record->d.id << endl;
                    exit(1);
                }
                hostGenotypes[hostSamples[i]].push_back(allele_index);
            }
        }

        // Loop over symbiont records (a,A=sym, b,B=host)
        float a = 0.0, A = 0.0, b = 0.0, B = 0.0, AB = 0.0, freq;
        map<string, int> symGenos;
        map<string, int>::iterator its;
        map<string, vector<int>>::iterator ith;
        for(size_t k=0; k<symbiontGenotypes.size(); k++) {
            // Loop over strains
            symGenos = symbiontGenotypes[k];
            for(const auto& strain: symbiontSamples) { 
                its = symGenos.find(strain);
                if(its == symGenos.end()) continue; // Missing genotype

                // Loop over hosts, assuming equal frequency of each host genotype
                for(const auto& host: hostSamples) {
                    ith = hostGenotypes.find(strain);
                    if(ith == hostGenotypes.end()) continue; // Missing genotype

                    freq = strainFrequencies[make_pair(strain, host)];

                    if(its->second == 0) {
                        a += freq;
                    } else {
                        A += freq;
                    }

                    for(const auto& allele: ith->second) {
                        if(allele == 0) {
                            b += 1;
                        } else {
                            B += 1;
                            if(its->second == 1) {
                                AB += freq;
                            }
                        }
                    }

                }

            }

            // Adjust allele frequencies for missing data
            float h = b + B;
            b /= h;
            B /= h;
            float s = a + A;
            a /= s;
            A /= s;
            AB /= s;

            // Calculate LD
            cout << symbiontRs[k] << "\t" << hostRs << "\t" <<
                calcR2(AB, A, B) << endl;

        }
        hostGenotypes.clear();

    }

    // Destory header object, record object, file object, and genotype
    // reading objects
    bcf_destroy(record);
    bcf_hdr_destroy(hostHdr);
    bcf_close(hostVCF);
    free(gt_arr);

}
