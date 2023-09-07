import * as fs from "fs";
import * as yaml from "js-yaml";
import { OpenAI } from "langchain/llms";
import { JsonSpec } from "langchain/tools";
import { JsonToolkit, createJsonAgent } from "langchain/agents";
import 'dotenv/config'

export const run = async () => {
  let data;
  try {
    const yamlFile = fs.readFileSync("data/case_ACTIVSg200 copy point.json", "utf8");
    data = yaml.load(yamlFile);
    if (!data) {
      throw new Error("Failed to load OpenAPI spec");
    }
  } catch (e) {
    console.error(e);
    return;
  }
  ////////////////// console.log(new JsonSpec(data))
  const toolkit = new JsonToolkit(new JsonSpec(data));
  
  const model = new OpenAI({modelName: "text-davinci-003", openAIApiKey: process.env.OPENAI_API_KEY, temperature: 0 });
  const executor = createJsonAgent(model, toolkit);

  const input = `What keys exist in the json file?`
  console.log(`Executing with input "${input}"...`);

  const result = await executor.call({ input });

  console.log("result output: \n")
  console.log(`Got output ${result.output}`);

  console.log("results: \n")
  console.log(result)

  console.log("intermediate: \n")
  console.log(
    `Got intermediate steps ${JSON.stringify(
      result.intermediateSteps,
      null,
      2
    )}`
  );
};

run(); 