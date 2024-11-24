# file name: UKBB_filtering
#
# year: 2023
# month: 10
# day: 09
#
library(susieR)
library(LDlinkR)
library(stringr)
#
#-------------------------------------------------------------------------------
# Statistical settings and input
#-------------------------------------------------------------------------------
#
pco_c<-5e-08 # p value cut-off for common variants
pco_uc<-1e-07 # p value cut-off for uncommon variants
pco_HWE<-1e-06 # p value for HW equilibrium
MAFco_c<-0.05 # lower cut-off for minor allele frequency of common variants
MAFco_uc<-0.01 # lower cut-off for minor allele frequency of uncommon variants
#
population<-seq(1:8) # 1 - 20002_1482 females
                      # 2 - 20002_1482 males
                      # 3 - 20002_1482 both
                      # 4 - R53 females
                      # 5 - R53 males
                      # 6 - R53 both
                      # 7 - 20002_1542 females
                      # 8 - 20002_1542 both
#
for (pop in population) {
  select<-pop
  if (select==1) {
    phenotype<-"20002_1482" # select phenotype
    sex<-"female" # select sex group
    temp<-read.csv("phenotypes.female.tsv",sep="\t")
    temp<-subset.data.frame(temp,phenotype=="20002_1482")
    n_cases<-temp$n_cases
    n_controls<-temp$n_controls
  }
  if (select==2) {
    phenotype<-"20002_1482" # select phenotype
    sex<-"male" # select sex group
    temp<-read.csv("phenotypes.male.tsv",sep="\t")
    temp<-subset.data.frame(temp,phenotype=="20002_1482")
    n_cases<-temp$n_cases
    n_controls<-temp$n_controls
  }
  if (select==3) {
    phenotype<-"20002_1482" # select phenotype
    sex<-"both_sexes" # select sex group
    temp<-read.csv("phenotypes.both_sexes.tsv",sep="\t")
    temp<-subset.data.frame(temp,phenotype=="20002_1482")
    n_cases<-temp$n_cases
    n_controls<-temp$n_controls
  }
  if (select==4) {
    phenotype<-"R53" # select phenotype
    sex<-"female" # select sex group
    temp<-read.csv("phenotypes.female.tsv",sep="\t")
    temp<-subset.data.frame(temp,phenotype=="R53")
    n_cases<-temp$n_cases
    n_controls<-temp$n_controls
  }
  if (select==5) {
    phenotype<-"R53" # select phenotype
    sex<-"male" # select sex group
    temp<-read.csv("phenotypes.male.tsv",sep="\t")
    temp<-subset.data.frame(temp,phenotype=="R53")
    n_cases<-temp$n_cases
    n_controls<-temp$n_controls
  }
  if (select==6) {
    phenotype<-"R53" # select phenotype
    sex<-"both_sexes" # select sex group
    temp<-read.csv("phenotypes.both_sexes.tsv",sep="\t")
    temp<-subset.data.frame(temp,phenotype=="R53")
    n_cases<-temp$n_cases
    n_controls<-temp$n_controls
  }
  if (select==7) {
    phenotype<-"20002_1542" # select phenotype
    sex<-"female" # select sex group
    temp<-read.csv("phenotypes.female.tsv",sep="\t")
    temp<-subset.data.frame(temp,phenotype=="20002_1542")
    n_cases<-temp$n_cases
    n_controls<-temp$n_controls
  }
  if (select==8) {
    phenotype<-"20002_1542" # select phenotype
    sex<-"both_sexes" # select sex group
    temp<-read.csv("phenotypes.both_sexes.tsv",sep="\t")
    temp<-subset.data.frame(temp,phenotype=="20002_1542")
    n_cases<-temp$n_cases
    n_controls<-temp$n_controls
  }
  if (select==9) {
    phenotype<-"I95" # select phenotype
    sex<-"female" # select sex group
    temp<-read.csv("phenotypes.female.tsv",sep="\t")
    temp<-subset.data.frame(temp,phenotype=="I95")
    n_cases<-temp$n_cases
    n_controls<-temp$n_controls
  }
  if (select==10) {
    phenotype<-"I95" # select phenotype
    sex<-"male" # select sex group
    temp<-read.csv("phenotypes.male.tsv",sep="\t")
    temp<-subset.data.frame(temp,phenotype=="I95")
    n_cases<-temp$n_cases
    n_controls<-temp$n_controls
  }
  if (select==11) {
    phenotype<-"I95" # select phenotype
    sex<-"both_sexes" # select sex group
    temp<-read.csv("phenotypes.both_sexes.tsv",sep="\t")
    temp<-subset.data.frame(temp,phenotype=="I95")
    n_cases<-temp$n_cases
    n_controls<-temp$n_controls
  }
  #
  input_phenotype<-gsub(" ","",paste(phenotype,".gwas.imputed_v3.",sex,".tsv.bgz"))
  #
  # Read the input file with the results for the phenotype/sex of choice
  #
  bgz<-gzfile(input_phenotype)
  mydata<-read.csv(bgz,header=T,sep="\t")
  #
  #-------------------------------------------------------------------------------
  # Search for variants significantly different between patients and controls
  #-------------------------------------------------------------------------------
  #
  # Filter common variants by p value and MAF, remove low confidence variants
  #
  v_ph_c<-subset.data.frame(mydata,minor_AF>=MAFco_c&pval<=pco_c&low_confidence_variant=="false")
  #
  # Filter uncommon variants by p value and MAF, remove low confidence variants
  #
  v_ph_uc<-subset.data.frame(mydata,minor_AF>=MAFco_uc&minor_AF<MAFco_c&
                               pval<=pco_uc&low_confidence_variant=="false")
  #
  remove(mydata) # free memory
  #
  # Read a version of file variants.tsv with only the columns "variant" and p_hwe
  #
  mydata<-read.csv("variants_p_hwe.tsv",header=T,sep="\t")
  #
  # Remove from v_ph_c the variants outside hwe
  #
  v_ph_c<-merge(v_ph_c,mydata,by="variant")
  #
  # Remove from v_ph_uc the variants outside hwe
  #
  v_ph_uc<-merge(v_ph_uc,mydata,by="variant")
  #
  remove(mydata) # free memory
  #
  # Calculate allele frequency
  #
  AAF_cases<-v_ph_c$ytx/(n_cases*2)
  AAF_controls<-(v_ph_c$AC-v_ph_c$ytx)/(n_controls*2)
  v_ph_c$Alt_AF_cases<-AAF_cases
  v_ph_c$Alt_AF_controls<-AAF_controls
  #
  AAF_cases<-v_ph_uc$ytx/(n_cases*2)
  AAF_controls<-(v_ph_uc$AC-v_ph_uc$ytx)/(n_controls*2)
  v_ph_uc$Alt_AF_cases<-AAF_cases
  v_ph_uc$Alt_AF_controls<-AAF_controls
  #
  #-------------------------------------------------------------------------------
  # Fine mapping
  #-------------------------------------------------------------------------------
  #
  # Function for the correlation matrix
  #
  Cor_Matrix<-function(n_v,f,Alt_AF) {
    errorD<-0
    D<-matrix(0,nrow=n_v,ncol=n_v)
    R<-matrix(1,nrow=n_v,ncol=n_v)
    for (j in 1:(n_v-1)) {
      if (j<n_v) {
        for (h in (j+1):n_v) {
          D[j,h]<- f[j,h] - Alt_AF[j]*Alt_AF[h]
          f_a<-min(Alt_AF[j], 1-Alt_AF[j])
          f_b<-min(Alt_AF[h], 1-Alt_AF[h])
          if (f_a<f_b) {                                     # test the range for D
            sup=(1-f_b)*f_a
            inf=-f_a*f_b
          }else { 
            sup=(1-f_a)*f_b
            inf=-f_a*f_b  
          }
          if (D[j,h]>sup) {
            D[j,h]<-sup
            errorD<-errorD+1
          }
          if (D[j,h]<inf) {
            D[j,h]<-inf
            errorD<-errorD+1
          }
          R[j,h]<-D[j,h]/sqrt( Alt_AF[j]*(1-Alt_AF[j])*Alt_AF[h]*(1-Alt_AF[h]) )
          if (R[j,h]>1) R[j,h]<-1                            # test the range for R
          if (R[j,h]<(-1)) R[j,h]<--1
        }
      }
    }
    for (j in 2:n_v) {
      for (h in 1:(j-1)) {
        R[j,h]<-R[h,j]
      }
    }
    Cor_Mat<-list()
    Cor_Mat[[1]]<-R
    Cor_Mat[[2]]<-errorD
    return(Cor_Mat)
  }
  #
  # Calculate the binary haplotype frequencies and the correlation matrix between
  # the variants to fine map
  #
  for (index in 1:2) {
    if (index==1) {
      df<-v_ph_c
      stringa<-".common"
    }
    if (index==2) {
      df<-v_ph_uc
      stringa<-".uncommon"
    }  
    n_v<-length(df[,1])
    if (n_v>1) {
      #
      # Interrogate LDhap for haplotype frequencies
      #
      haplo<-LDhap(snps=df$rsid,pop="GBR",token="100e14bc74eb",genome_build="grch37",table_type="haplotype")
      var<-LDhap(snps=df$rsid,pop="GBR",token="100e14bc74eb",genome_build="grch37",table_type="variant")
      #
      Alt<-c()
      for (i in 1:n_v) Alt[i]<-str_extract(gsub(" ","",var[i,3]),"(?<=,)[^,]") # alternative alleles
      #
      tot<-sum(as.numeric(haplo$Count))
      haplo$Frequency<-as.numeric(haplo$Count)/tot # more precise calculation of Frequency
      n_h<-length(haplo[,1])
      #
      f<-matrix(0,nrow=n_v,ncol=n_v)
      for (j in 1:n_v) {
        for (i in 1:n_h) {
          if (j<n_v) {
            for (h in (j+1):n_v) {
              if (haplo[i,j]==Alt[j]) {
                if (haplo[i,h]==Alt[h]) {
                  f[j,h]<-f[j,h]+haplo$Frequency[i]
                }
              }
            }
          }
        }
      }
      #
      # Calculate the correlation matrix with allele frequency from GBR population
      # and alternative allele frequencies from the control group
      #
      Cor_Mat<-Cor_Matrix(n_v,f,df$Alt_AF_controls)
      #
      # PIP calculation with SusieR
      #
      z_scores = df$beta/df$se # we calculate the z score
      #
      R = data.matrix(Cor_Mat[[1]])
      #
      # Set the max number admitted of causal variants (L) and the sample size (n) 
      #
      n<-n_cases+n_controls # sample size
      phi<-n_cases/n_controls
      ne<-n*phi*(1-phi) # effective sample size
      L<-1
      #
      # Use the function fitted_rss of the library susieR to calculate 
      # the posterior inclusion probability (PIP) of each variant
      #
      fitted_rss = susie_rss(z_scores, R, n = ne, L = L, z_ld_weight = 0)
      #
      # Plot PIP as a function of the variant (numbered from 1 to length(beta))
      #
      susie_plot(fitted_rss, y="PIP",xlab="variants")
      #
      # Write on the screen a summary of the results
      #
      SusieResult<-summary(fitted_rss)
      #
      PIP<-SusieResult[[1]][order(SusieResult[[1]]$variable),]
      #
      df$CS<-PIP$cs
      df$PIP<-PIP$variable_prob
      #
      # Joint model approximation
      #
      sc<-sqrt(2*df$minor_AF*(1-df$minor_AF)) # scaling factor
      beta.sc<-df$beta*sc # scaled regression coefficients for marginal models
      lambda.sc<-solve(R,beta.sc) # scaled regression coefficients for joint model
      lambda<-lambda.sc/sc # regression coefficients for joint model on allelic scale
      SE.lambda.sc<-sqrt(diag(solve(R))/ne) # scaled standard errors for joint model
      SE.lambda<-SE.lambda.sc/sc # standard errors for joint model on allelic scale
      z_score.lambda<-lambda/SE.lambda
      pval.lambda<-2*pnorm(-abs(z_score.lambda),0,1,lower.tail=T)
      #
      df$lambda<-lambda
      df$se_lambda<-SE.lambda
      df$z_lambda<-z_score.lambda
      #
    }
    #
    # write the output of filtering and fine mapping
    #
    file_name<-gsub(" ","",paste(phenotype,".",sex,stringa,".tsv"))
    write.csv(df,file=file_name)
  }
}


