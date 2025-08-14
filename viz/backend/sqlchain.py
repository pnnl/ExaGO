from langchain_openai import ChatOpenAI
from langchain_community.utilities import SQLDatabase
from langchain_experimental.sql import SQLDatabaseChain
import config
from langchain.output_parsers import CommaSeparatedListOutputParser
from langchain.prompts import PromptTemplate, ChatPromptTemplate, HumanMessagePromptTemplate
import sqlalchemy as sqldb
from sqlalchemy import text


def sqlchain(input_text):
    llm = ChatOpenAI(openai_api_key=config.openai_key,
                 model_name="gpt-4", temperature=0, verbose=True)
    db = SQLDatabase.from_uri(
        f"postgresql+psycopg2://postgres:{config.sql_key}@localhost:5432/{config.database_name}")
    mydb = sqldb.create_engine(
        f"postgresql+psycopg2://postgres:{config.sql_key}@localhost:5432/{config.database_name}")
    myconnection = mydb.connect()

    _CUSTOMIZE__TEMPLATE = """You are a PostgreSQL expert. Given an input question, first create a syntactically correct PostgreSQL query to run then look at the results of the query and return the answer to the input question.
    You must always query for the name column (e.g., generation_name, line_name, bus_name). Never query for all columns from a table.  Wrap each column name in double quotes (") to denote them as delimited identifiers.

    When users query about state or county, you can use ST_GeomFromText(wkt) and postgis functions to calculate spatial relationship between geo entities. 
    When generating SQL queries involving spatial data stored for states and counties in WKT format (e.g., in columns like coordinates or wkt), always cast the WKT string to a PostGIS geometry using ST_GeomFromText(column, 4326).
    For example, ST_GeomFromText(wkt, 4326)

    If checking if a point lies within a polygon (e.g., a generator in a state or county), use:
    ST_Contains(ST_GeomFromText(state.wkt, 4326), ST_GeomFromText(generator.coordinates, 4326))


    Pay attention to use only the column names you can see in the tables below. Be careful to not query for columns that do not exist. Also, pay attention to which column is in which table.
    First look the postgre database for an answer, if you can't find related answer from the database, you can use you own knowlwedge to answer the questions.

    Use the following format:

    Question: Question here
    SQLQuery: SQL Query to run
    SQLResult: Result of the SQLQuery
    Answer: text answer

    For the text answer,
      - Use appropriate formatting that shows the answer clearly.
      - Include a short summary for the text answer for more than 5 data points.
      - Use appropriate units (for e.g. MW, KV, or pu) as needed.

    """
    MY_PROMPT_SUFFIX = """Only use the following tables:
    {table_info}

    Question: {input}"""

    MY_POSTGRES_PROMPT = PromptTemplate(
        input_variables=["input", "table_info"],
        template=_CUSTOMIZE__TEMPLATE + MY_PROMPT_SUFFIX,
    )

    text_chain = SQLDatabaseChain.from_llm(
        llm, db, verbose=True, return_intermediate_steps=True, prompt=MY_POSTGRES_PROMPT)
    # sql_chain = SQLDatabaseChain.from_llm(llm, db, verbose=True,return_intermediate_steps=True, return_direct=True)

    # ResultSet = tempr.fetchall()

    # format input and output with prompt template
    # output_parser = CommaSeparatedListOutputParser()

    # format_instructions = output_parser.get_format_instructions()

    # _DEFAULT_TEMPLATE = """Given an input question, first create a syntactically correct query to run. The return of the query should always include the 'id' field. Then look at the results of the query and return the answer.

    # # Use the following format:

    # # Question: "Question here"
    # # SQLQuery: "SQL Query to run"
    # # SQLResult: "Result of the SQLQuery"
    # # Answer: ""

    # Question: {input}"""
    # prompt = PromptTemplate(
    #     input_variables=["input"], template=_DEFAULT_TEMPLATE
    # )

    # prompt = PromptTemplate(
    # template="Answering the following questions{query}.\n{format_instructions}",
    # input_variables=["query"],
    # partial_variables={"format_instructions": format_instructions}
    # )
    # input_query = "How many generations in total are there in Illinois?"
    # _input = prompt.format(input=input_text)
    # "show me all the wind generations with capacity higher than 100 in Illinois"
    # Find the top 3 closest generations to the generation 'SPRINGFIELD 5' with capacity higher than 100 in Illinois.

    text_result = ""
    query_dict = []
    try:
        output = text_chain.invoke(input_text)
        # sql_result = output["intermediate_steps"][3]
        text_result = output['result']
        # print(output["intermediate_steps"])
        sql_cmd = output["intermediate_steps"][1]

        # list_output = output_parser.parse(output)

        tempr = myconnection.execute(text(sql_cmd)).fetchall()

        query_dict = [dict(record._mapping) for record in tempr]
    except Exception as error:
        error_message = str(error)
        if "maximum" in error_message and "length" in error_message:
            text_result = "Check the visualization for updated results"
            sql_cmd = error.intermediate_steps[1]
            tempr = myconnection.execute(text(sql_cmd)).fetchall()
            query_dict = [dict(record._mapping) for record in tempr]
        elif "Rate" in error_message and "limit" in error_message:
            text_result = "Check the visualization for updated results"
            sql_cmd = error.intermediate_steps[1]
            tempr = myconnection.execute(text(sql_cmd)).fetchall()
            query_dict = [dict(record._mapping) for record in tempr]
        else:
            print("error:")
            print(error)
            text_result = "Sorry, I can't find the answer to your question"
            query_dict = []

    # print(query_dict)
    print(text_result)
    print(sql_cmd)
    return {
        "text": text_result,
        "result_list": query_dict
    }
# print(list_output)

# sqlchain('show me US generations that are in Texas state?')
