# Import flask and datetime module for showing date and time
from flask import Flask
import datetime
from flask_cors import CORS
from flask.globals import request
from sqlchain import sqlchain
from flask import jsonify


x = datetime.datetime.now()

# Initializing flask app
app = Flask(__name__)
CORS(app)


# Route for seeing a data
@app.route('/data', methods=['POST'])
def get_time():
    input_text = request.get_json(force=True)["inputText"]
    # print(request.json)
    print(input_text)
    # outputtext = sqlagent(input_text)
    output = sqlchain(input_text)
    # Returning an api for showing in  reactjs
    return jsonify(output)


# Running app
if __name__ == '__main__':
    app.run(debug=True)
