# tamipami:  a tool for finding the TAM and PAM sites of novel endonucleases

This tools takes a control fastq sequencing library and an experimental library treated with 
the endonuclease you are trying to find a TAM /PAM site for.  The two libraries are generated 
from a pool of sequences containing a guide site and a mixture of random bases flanking the 
guide region. When a endonuclease finds the PAM/TAM site and a guide RNA it cuts and those sequences are 
depleted in the sample library relative to the control.


# Web Application 

To launch a the streamlit web application locally at http://localhost:8501 run this commmand:

```{bash}
streamlit run app.py 
```
