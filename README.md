)# Modeling-Gender-Wise-Gene-Response-to-Smoking

**Introduction **

The motivation on is to identify genes that exhibit differential response to smoke based on gender. A two-way analysis of variance (ANOVA) framework is utilized to generate p-values for each row in the dataset. Thep-values are then analysed by drawing a histogram to gain insights into the significance of the findings.

**Implementation **

In the starting ,the necessary libraries, such as pandas, numpy, re, scipy.stats, and matplotlib.pyplot are imported. The dataset is read using pandas from the 'Raw Data_GeneSpring.txt' file, assuming it is located in the 'data' folder relative to the Python file. A regular expression pattern is defined to identify the columns in the dataset that match the specified format. This is done so that only the ‘GSM’ genes are captured and the gender and smoking status is being group filtered to make the data corresponding to gender and smoking status. The code then iterates over the columns to determine the count of columns 
with the desired format. 

Next, the code initializes matrices and dictionaries required for computations. I utilise these matrices for calcula on of f-sta s c and ul mately p-values for each row in the data. 
For each row, the code extracts the relevant data and populates matrices for gender(Male/Female) and smoking status (Smoker/Non-smoker). The rank of matrices D and N is determined, and the necessary calculations are performed to obtain the F-statistic, p-value, and degrees of freedom. 

The computed p-values are stored in a dictionary stat_dict along with the row index and F-statistic. Additionally, a count of rows with p-values less than or equal to 0.05 is obtained. Finally, the code defines a plo ng func on and calls it to generate histograms of the p-values using different numbers of bins.

**Histogram of p-values **

To visualize the distribu on of p-values, histograms were generated with varying numbers of  bins (10, 20, 100, 1000, and the total number of rows~41k).The varying number of bins help observe the p-values in a unique way.  The histograms provide insights into the statistical significance of the findings and the overall distribution of p-values. 

**Results**

![image](https://github.com/Aksheit-Saxena/Modeling-Gender-Wise-Gene-Response-to-Smoking/assets/58588004/5d4cb2dc-bc4a-44fd-8579-81847a545572)


![image](https://github.com/Aksheit-Saxena/Modeling-Gender-Wise-Gene-Response-to-Smoking/assets/58588004/d922cf26-23ae-4ab4-bee4-d2ef0a04b434)



![image](https://github.com/Aksheit-Saxena/Modeling-Gender-Wise-Gene-Response-to-Smoking/assets/58588004/e2ff2da1-2e98-4dca-b55f-7bb21bd9085e)


**Conclusion**

The work successfully u lized a two-way ANOVA framework to iden fy genes that respond differently to smoke based on gender. There are six different histogram plo ed for varied number of bins to capture how the p-values are spread. In the initial few graphs it is not clear which p-values take preference but in the later graphs it is clear that the p-value 
around 0.3 has the highest occurrence. By applying a threshold of 0.05, the number of rows showing significant differences in gene response to smoke in men and women was 
determined by the variable ‘fiveper_pvalue’  in code which has 811 rows.



