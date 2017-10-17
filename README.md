# ADHD_classification_DBN_extraction
Guide : Prof Sundaram Suresh (NTU- Singapore) Area: Deep learning neural networks for feature extraction in high dimensional neuro imaging data. Tools used: Standard neuro imaging software for preprocessing, a MATLAB deep learning toolbox DeeBNet. I used deep learning algorithms including RBM’s and CNN’s to train on an open source MRI data set and classify unseen fMRI scans as having ADHD or not. I was able to achieve accuracy scores of 64% which is incrementally better than the current start of art(as of 2016). The project was challenging due to the high dimensionality of the input data and the meager number of test samples.

Since the clinical fmri data is high-dimensional naive, shallow classifiers are not able to do a good job. 
The basic idea behind the project was to investigate if RBM's/ DBN's can extract useful higher level features.

I implemented a RBM based multilayer Deep Belief Network and extracted features from the raw dataset. I had to 
perform a lot of experiments to train the network properly. 

The results indicate that the features extracted using the above lead to significantly better classification results.
I was able to achieve results as high as 75% compared to the previous baseline of 65%.

The .docx files are the reports which have more details. 
