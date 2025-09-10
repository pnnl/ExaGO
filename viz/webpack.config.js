// NOTE: To use this example standalone (e.g. outside of deck.gl repo)
// delete the local development overrides at the bottom of this file

const Dotenv = require('dotenv-webpack');

const CONFIG = {
  mode: 'development',

  entry: {
    app: './app.js'
  },

  output: {
    library: 'App',
    publicPath: '/',
    // sourceMapFilename: '[name].js.map'
  },

  devServer: {
    historyApiFallback: {
      disableDotRule: true,
      index: '/index.html'
    },
    contentBase: __dirname,
    compress: true,
    port: 8080,
    watchOptions: {
      poll: true
    },
    publicPath: '/',
    hot: true
  },

  devtool: 'source-map',

  plugins: [
    new Dotenv()
  ],

  module: {
    rules: [
      {
        // Transpile ES6 to ES5 with babel
        // Remove if your app does not use JSX or you don't need to support old browsers
        test: /\.js$|jsx/,
        loader: 'babel-loader',
        exclude: [/node_modules\/(?!(firebase|@firebase|es-toolkit)\/).*/],
        options: {
          presets: [
            ['@babel/preset-env', {
              targets: {
                browsers: ['last 2 versions', 'ie >= 11']
              }
            }],
            '@babel/preset-react'
          ],
          plugins: [
            '@babel/plugin-proposal-optional-chaining',
            '@babel/plugin-proposal-nullish-coalescing-operator'
          ]
        }
      }
    ]
  },
  resolve: {
    extensions: ['.js', '.jsx']
  }

};

// Export the configuration directly
module.exports = CONFIG;
