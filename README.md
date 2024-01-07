# Decide
a web visualizer to display chemical compounds and associated data

1. Clone git

```git clone https://github.com/maximestephan/decide.git```

2. Create environment

on Windows:
    ```py -m venv decideenv```
on Linux/macOS:
    ```python3 -m venv decideenv```

3. Activate environment

on Windows:
    ```.\decideenv\Scripts\activate```
on Linux/macOS: 
    ```./decideenv/bin/activate```
    
( if you get a permission deny, execute this command 
    ```chmod u+x ./decideenv/bin/activate```
)

4. Install requirements

```pip install -r requirements.txt```

5. run flask

```flask --app main.py run --debug -p 8080```

6. run Decide locally on your browser

open a browser on http://127.0.0.1:8080