from langchain import OpenAI, SQLDatabase, SQLDatabaseChain
import config
from langchain.output_parsers import CommaSeparatedListOutputParser
from langchain.prompts import PromptTemplate, ChatPromptTemplate, HumanMessagePromptTemplate
import sqlalchemy as sqldb
from sqlalchemy import text


def sqlchain(input_text):
    llm = OpenAI(openai_api_key=config.openai_key,
                 model_name="gpt-3.5-turbo", temperature=0, verbose=True)
    db = SQLDatabase.from_uri(
        f"postgresql+psycopg2://postgres:{config.sql_key}@localhost:5432/US west power grid")
    mydb = sqldb.create_engine(
        f"postgresql+psycopg2://postgres:{config.sql_key}@localhost:5432/US west power grid")
    myconnection = mydb.connect()

    _CUSTOMIZE__TEMPLATE = """You are a PostgreSQL expert. Given an input question, first create a syntactically correct PostgreSQL query to run then look at the results of the query and return the answer to the input question.
    You must always query for the name column (e.g., generation name, line name, bus name). Never query for all columns from a table.  Wrap each column name in double quotes (") to denote them as delimited identifiers.
    When users query about state or county, you can use ST_GeomFromText(wkt) and postgis functions to calculate spatial relationship between geo entities. 
    Pay attention to use only the column names you can see in the tables below. Be careful to not query for columns that do not exist. Also, pay attention to which column is in which table.
    First look the postgre database for an answer, if you can't find related answer from the database, you can use you own knowlwedge to answer the questions. 

    Use the following format:

    Question: Question here
    SQLQuery: SQL Query to run
    SQLResult: Result of the SQLQuery
    Answer: text answer

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
        output = text_chain(input_text)
        # sql_result = output["intermediate_steps"][3]
        text_result = output['result']
        # print(output["intermediate_steps"])
        sql_cmd = output["intermediate_steps"][1]

        # list_output = output_parser.parse(output)

        tempr = myconnection.execute(text(sql_cmd)).fetchall()

        query_dict = [dict(record._mapping) for record in tempr]
    except Exception as error:
        if (("maximum" in error.user_message) and ("length" in error.user_message)):
            text_result = "Check the visualization for updated results"
            sql_cmd = error.intermediate_steps[1]
            tempr = myconnection.execute(text(sql_cmd)).fetchall()
            query_dict = [dict(record._mapping) for record in tempr]
        elif (("Rate" in error.user_message) and ("limit" in error.user_message)):
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
