import { OpenAI } from "langchain/llms/openai";
import { SqlDatabase } from "langchain/sql_db";
import { createSqlAgent, SqlToolkit } from "langchain/agents";
import { DataSource } from "typeorm";
import 'dotenv/config'
// import { PromptTemplate } from "langchain/prompts";


/** This example uses Chinook database, which is a sample database available for SQL Server, Oracle, MySQL, etc.
 * To set it up follow the instructions on https://database.guide/2-sample-databases-sqlite/, placing the .db file
 * in the examples folder.
 */
 async function chatModel(inputText) {
    const datasource = new DataSource({
        type: "postgres",
        host: "localhost",
        port: 5432,
        username: "postgres",
        password: process.env.SQL_KEY,
        database: "postgres",
    });
    const db = await SqlDatabase.fromDataSourceParams({
        appDataSource: datasource,
    });
    //   console.log(db)
    const model = new OpenAI({ modelName: "gpt-3.5-turbo", openAIApiKey: process.env.OPENAI_API_KEY, temperature: 0 });
    const toolkit = new SqlToolkit(db, model);
    const executor = createSqlAgent(model, toolkit);

      const input =  inputText
    //   "How many generations are there in the table point? "
    //   `How many generations are there with pcap higher than 200 in the table "point"?`
    //   `Find the top 3 closest points to the point id 1 with "pcap" higher than 100 in the table "point". The point coordinate infomration is in the column wkt. You can use postgis functions.`
    //   `List id of generations with connection higher than 10000 in the table power1.`;

      console.log(`Executing with input "${input}"...`);

    // console.log('\n')
    // console.log(inputText)
    // const template = "{inputText}";
    // const prompt = new PromptTemplate({
    //     template: template,
    //     inputVariables: ["inputText"],
    // });
    // const promptText = await prompt.format({ inputText: inputText});
    // console.log(promptText)
    const result = await executor.call({ input });



    console.log(`Got output ${result.output}`);
    console.log('\n')

    console.log(
        `Got intermediate steps ${JSON.stringify(
            result.intermediateSteps,
            null,
            2
        )}`
    );

    await datasource.destroy();
    return result.output
};

// export {chatModel};
const inputText =  `List id and capacity of the top 5 generations that have the most capacity in table point200`
// "Return the most robust point in table point200?"
// `how many points have capacity higher than 200 in the table point200?`
//  `How many points do not have pcap higher than 200 in the table "point"?`
// `Point in which color generate the most power?`
// `Which point in green has the most pcap in table poin200?`
// `How many points in total in table point200? `
// `What is the closest point to point id 1 (except itself) in Eclidean distance in the table point200? The point coordinate infomration is in the column wkt. You can use postgis functions`
// `What are the names the id of generations are there with pcap higher than 200 in the table "point"?`
//   "How many points are there in the table point200? "
    //   `How many generations are there with pcap higher than 200 in the table "point"?`
    //   `Find the top 3 closest points to the point id 1 with "pcap" higher than 100 in the table "point". The point coordinate infomration is in the column wkt. You can use postgis functions.`
    //   `List id of generations with connection higher than 10000 in the table power1.`;
chatModel(inputText) 


