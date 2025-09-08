# Import flask and datetime module for showing date and time
from flask import Flask
import datetime
from flask_cors import CORS
from flask.globals import request
from flask import jsonify

# Temporarily disable database-dependent imports to get basic server running
# from sqlchain import sqlchain
# from agent import run_agent_json


x = datetime.datetime.now()

# Initializing flask app
app = Flask(__name__)
CORS(app)


# Route for basic health check
@app.route('/health', methods=['GET'])
def health_check():
    return jsonify({"status": "healthy", "timestamp": str(datetime.datetime.now())})


# Route for seeing a data
@app.route('/data', methods=['POST'])
def get_time():
    input_text = request.get_json(force=True)["inputText"]
    print(f"Received input: {input_text}")
    
    # Temporary response while database functionality is being set up
    output = {
        "response": f"Backend received your message: '{input_text}'. Database functionality is currently being configured.",
        "timestamp": str(datetime.datetime.now())
    }
    
    # TODO: Re-enable when database issues are resolved
    # output = run_agent_json(input_text)
    
    return jsonify(output)


# Running app
if __name__ == '__main__':
    # Run on all interfaces for Docker, port 5000
    app.run(host='0.0.0.0', port=5000, debug=True)
