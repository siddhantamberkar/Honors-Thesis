# SMILES-decoder
This project is for my UC Berkeley Department of Letters & Sciences Physics Honors Thesis. The goal is to create a framework to identify amines from SMILES strings with the purpose of using it towards data mining large chemical structure databases such as CSD. Amines are key in carbon capture technology and this work would assist in a blind crystal search of amines to match thermodnymical properties with the SMILES strings and ultimately identify which amines are key to efficient carbon capture. 

I first looked at polarizability of the molecules of interest (amines) since it lends key insight into the melting point of respective crystal structure. In carbon capture, we are interested in crystal strcutures with high melting points since they are capable of resisting material loss during adsorption/desorption cycles. 

To successfully run these calculations, I used the [ORCA](https://www.faccts.de/orca/) open-source quantum chemistry package. You'll also need to download ORCA locally and set the ORCA_PATH environment variable or add ORCA to your system PATH.
