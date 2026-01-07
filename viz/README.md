# ExaGO visualization (experimental)
ExaGO has an experimental visualization platform for visualizing the results of OPFLOW on a map provided the geospatial coordinates for the network are available. On launching the visualization, a webpage displays the given power system network overlayed on a geospatial map. The (experimental) visualization features include: 
  - Map-based network layout of the grid
  - Reading in the grid data through a geojson file.
  - Fly-in (zoom in) on bus, branch, county 
  - County-based load and voltage heatmaps
  - Bar chart layer and double pie-chart for generation dispatch and capacity
  - Filters based on network voltage, generation dispatch level, and voltage/load level
  - Zoomed-in display of county load or voltage (aggregated)


## Installation
ExaGO visualization uses the following tools to generate the visuals.
- [Node.js@v24.10.0](https://nodejs.org/es/blog/release/v24.10.0)
- Facebook's [React](https://reactjs.org/) framework
- Uber's [Deck.gl](https://deck.gl/docs) visualization
- [React-map-gl](https://visgl.github.io/react-map-gl/) framework
- [Chart.js](https://www.chartjs.org/)
- Yarn 1.22.22

Before launching the visualization, one needs to install these packages. This can be done with the following steps:

1. Install Node Version Manager (NVM). On MAC use `brew install nvm`.
2. Install [Node.js] version 24 (https://nodejs.org/en/) using `nvm install 24`
3. Select Node 24 by `nvm use 24`
4. Install Yarn. Do `npm install --global yarn`
5. Do `yarn install` in this directory (`viz`) to install all the dependencies.
6. Go to the `viz/backend` subdirectory and use the `pip install -r requirements.txt` command to install all the Python dependencies.


## Preparing input data files for visualization
The visualization uses a `JSON` formatted file as an input. This `JSON` file has a specific structure (To do: explain structure for the file) and there are several sample files for different network in the `data` subdirectory.
This input JSON file can be either created externally OR generated as an output of the `OPFLOW` application. When using OPFLOW, the following command will generate the input JSON file. The generated file will be name as `opflowout.json`.

```
./opflow -netfile <netfile> -save_output -opflow_output_format JSON -gicfile <gicfilename>  
```

Note that the `OPFLOW` application is available in the `$EXAGO_INSTALL/bin` directory where `$EXAGO_INSTALL` is the ExaGO installation directory.

The above command will run a `OPFLOW` on the given network and generate an output file called `opflowout.json`. The `-gicfile` is an additional option one can provide to provide the file that has the geospatial coordinates (latitude/longitude) for the network. If the geospatial coordinates are not provided then OPFLOW draws the network as a circle. It is highly recommended that one provides the geospatial coordinate file as an input to display the network correctly on the map. The geospatial coordinate file should have the same format as used for the [Electric Grid Test Case Repository](https://electricgrids.engr.tamu.edu/) synthetic networks. 

For example, with Texas 2000 bus synthetic data, executing the following `opflow` will produce the `opflowout.json` output. The case files are provided in the data folder.

```
opflow -netfile case_ACTIVSg2000.m -save_output -opflow_output_format JSON -gicfile ACTIVSg2000_GIC_data.gic
```

Next, go to the `viz` folder and run the following python script `geninputfile.py` from the `viz` folder to load the JSON file (`path/to/opflowout/opflowout.json`) in the visualization script. It will copy the `json` file to the `viz/data` subdirectory and create/overwrite a file named `viz/src/module_casedata.js`. The `module_casedata.js` file is an application source file to load the data file `opflowout.json`. Note, the visualization tool expects the file (`opflowout.json`) to be present in `viz/data` forlder, so it is copied by this script.

```
python geninputfile.py path/to/opflowout/opflowout.json
```

Now you are ready to launch the visualization now. 

Note: If you have already created or have the JSON file externally without running the `opflow` command as instructed above, simply run the `geninputfile.py` script using the above command.

## Launch visualization
To launch the visualization, run

```
yarn start
```

This will open a webpage (e.g. `http://localhost:5173/`) on the default browser with the visualization of the given network. If the network is large, it may take a while to load the visualization. The browser may show option to terminate or Wait and you should click on Wait button.


The figures show the visualization of the synthetic electric grid. The data for developing this visualization was created by merging the synthetic dataset for the [Eastern](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/activsg70k/), [Western](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/activsg10k/), and [Texas](https://electricgrids.engr.tamu.edu/electric-grid-test-cases/activsg2000/) interconnects from the [Electric Grid Test Case Repository](https://electricgrids.engr.tamu.edu/)
### 2D synthetic US western grid network display
![](images/network_viz.PNG)

### 2D synthetic US western grid tranmission line flow display
![](images/flow_viz.PNG)

### 2.5D synthetic US western grid network display with generation overlapped and doughnut chart for generataion mix
![](images/generation_viz.PNG)

### 2.5D synthetic US western grid displaying load profile by counties
![](images/load_viz.PNG)

### 2.5D synthetic US western grid displaying network, flow, generation, and load
![](images/all_viz.PNG)

### Demo 
See [Here](../tutorials/demo1.ipynb) 

## ChatGrid
ChatGrid is a natural language query tool for ExaGO visualizations. It is powered by OpenAI GPT-3.5-Turbo and Langchain. ChatGrid allows users to query on ExaGO visualizations through natural language and returns text summaries and visual outputs as answers. The following flow chart shows the architecture design of ChatGrid.
![](images/chatgrid_arch.png)

### Dependencies
ChatGrid is built upon the following services and tools. 

- [OpenAI LLMs](https://platform.openai.com/docs/models/overview)
- [Langchain@0.0.233](https://python.langchain.com/docs/get_started/introduction.html) framework
- [PostGreSQL](https://www.postgresql.org/download/) database
- [Flask](https://flask.palletsprojects.com/en/2.3.x/) framework 

### Preparing the database for query
Behind the scenes, LLM translates natural language queries into SQL queries to retrieve data from a database. As the power grid network is a typical geospatial dataset, we choose PostgreSQl + PostGIS database for the convenience of conducting spatial queries. Please follow the steps below to set up your PostgreSQL database that contains the power grid network dataset. 

1. Convert data formats. 

    First, we need to convert the ExaGO output `.json` files to `.csv` files. The difference between the two data formats is that JSON stores attributes and values as dictionary pairs but CSV stores attributes and values as tables. You can write your own script for this conversion or use the provided script. 
    
      * Go to the `viz/backend` subdirectory and use the `pip install -r requirements.txt` command to install all the Python dependencies if already not done in previous steps. (Note: These steps are tested with Python 3.13.)


    To use the provided script, ensure that the ExaGO output `.json` file is in the `viz/data` subdirectory (if not already performed in the previous steps) and run the following script in the `viz/backend` subdirectory (replace the example filename with your json filename). This will create three CSV files: `generation.csv`, `bus.csv`, and `tranmission_line.csv`. We are assuming `opflowout.json` is the data json file present in `viz/data` folder.
            
    ```
    cd viz/backend
    python ../data/jsontocsv.py ../data/opflowout.json
    ```
    
    Now there should be 5 CSV files in `viz/backend` folder: 
        
      * `bus.csv`               
      * `generation.csv`        
      * `us_states.csv`
      * `counties.csv`          
      * `transmission_line.csv`
    
    
2. Download PostgreSQL database from this [link](https://www.postgresql.org/download/) and install it.

    * For MAC using brew you can install postgresql 14 using: `brew install postgresql@14`
    * Start the postgressql service: `brew services start postgresql@14`
    * Create a role: `psql -U "$USER" -d postgres`
    * If your username (`$USER`) starts with numbers `psql` will show error. In this case replace `$USER` to some other username. 
    * Execute the create role query: `CREATE ROLE postgres WITH LOGIN SUPERUSER PASSWORD 'ExaGO.2025';` Here `ExaGO.2025` is a password. Change to your preference.
      * If you forget the password, you can update it by running `ALTER USER user_name WITH PASSWORD 'new_password';` (replace `user_name` and `new_password`)
    * Exit to shell by entering `quit` and hitting Enter.
    * From command prompt type: `psql -U postgres -d postgres` If it works and you are in `psql` shell you are done. Exit from the shell using `quit`.


3. Create a PostgreSQL database and import the `.csv` files to it.

    There are many ways to do that. If you are not familiar with SQL commands, you can use `pgAdmin 4` (a GUI management program of PostgreSQL database, usually installed along with PostgreSQL) and follow [this instruction](https://learnsql.com/blog/how-to-import-csv-to-postgresql/) to import the `.csv` files (e.g., generation.csv, bus.csv, tranmission_line.csv, counties.csv, US states.csv) to PostgreSQL database.

    Note: 

      a. Please write down your username, password, database name and keep them in a safe place. You will need this information to connect to the database later. 

      b. Please be informative and accurate about your table names, and attribute names. Because this information can help LLM understand the dataset and performs better when dealing with user queries.

      c. Include US state and county information in your database to support spatial queries that related to state or county.
      
      d. To enter the CSV files into database using command prompt do: `PGPASSWORD=ExaGO.2025 ./create_db.sh  --db exago_db --schema-sql ./schema.sql --drop  --truncate`. Here `exago_db` is the database name. Use it in the configuration `config.py` file.
      
        * If you used different user name than `$USER` then use the following command: `PGUSER=user_name PGPASSWORD=ExaGO.2025 ./create_db.sh  --db exago_db --schema-sql ./schema.sql --drop  --truncate` (replace `user_name` with your user name.)
      
      e. This will create a database named `exago_db` with password `ExaGO.2025`. This information will be used to update the `config.py` file.
    

4. Connect to your database.

    Open the `config.py` file in the `viz/backend` subdirectory and replace `YOUR_DATABASE_PASSWORD` and `YOUR_DATABASE_NAME`  with your own database password and database name.


### Getting your OpenAI API key
ChatGrid uses GPT models from OpenAI to process natural language queries. To use LLMs from OpenAI, you first need to go to [OpenAI's Platform website](https://platform.openai.com) and sign in with an OpenAI account. Click your profile icon at the top-right corner of the page and select "View API Keys." Click "Create New Secret Key" to generate a new API key. 

Open the `config.py` in the `viz/backend` subdirectory replace `YOUR_OPENAI_KEY` with your own OpenAI API key.


<!-- data script -->
### Launch backend
ChatGrid uses Flask to host the service of receiving user queries and returning the data output and text summaries to update the visualizations on the frontend. Please follow the steps below to run the backend server.


* Run the following command in the `viz/backend` subdirectory
    ```
    python server.py
    ```
    
    This will start the backend server for receiving user queries, processing it by LLM and returning data outputs to the frontend. 

Now open the chat window on the frontend, type your queries, and enjoy ChatGrid!
![](images/chatgrid_case.png)

### Optional: model configuration 
If you would like to test different LLMs with ChatGrid, you can specify the `model_name="YOUR_LLM_MODEL"` in the `viz/backend/sqlchain.py` file. 
