import React, { useState } from 'react';
import { useAuth } from './AuthContext';
import { Box, Button, TextField, Typography, Paper, Alert, Stack, ThemeProvider } from '@mui/material';
import PersonIcon from '@mui/icons-material/Person';
import LockIcon from '@mui/icons-material/Lock';
import { authTheme } from './theme';

const Login = () => {
    const [email, setEmail] = useState('');
    const [password, setPassword] = useState('');
    const [confirmPassword, setConfirmPassword] = useState('');
    const [displayName, setDisplayName] = useState('');
    const [isSignUp, setIsSignUp] = useState(false);
    const [error, setError] = useState('');
    const [loading, setLoading] = useState(false);
    const { login, signup } = useAuth();

    async function handleSubmit(e) {
        e.preventDefault();

        if (isSignUp && password !== confirmPassword) {
            return setError('Passwords do not match');
        }

        try {
            setError('');
            setLoading(true);

            if (isSignUp) {
                await signup(email, password, displayName);
            } else {
                await login(email, password);
            }
        } catch (error) {
            setError(isSignUp ? 'Failed to create account: ' + error.message : 'Failed to log in: ' + error.message);
        }

        setLoading(false);
    }

    const toggleMode = () => {
        setIsSignUp(!isSignUp);
        setError('');
        setPassword('');
        setConfirmPassword('');
        setDisplayName('');
    };

    return (
        <ThemeProvider theme={authTheme}>
            <Box
                sx={{
                    display: 'flex',
                    justifyContent: 'center',
                    alignItems: 'center',
                    minHeight: '100vh',
                    background: 'linear-gradient(135deg, #74a9cf 0%, #a8c8e1 100%)',
                    padding: 2
                }}
            >
                <Paper
                    elevation={10}
                    sx={{
                        padding: 4,
                        borderRadius: 3,
                        width: '100%',
                        maxWidth: 450,
                        backgroundColor: 'rgba(255, 255, 255, 0.95)',
                        backdropFilter: 'blur(10px)'
                    }}
                >
                    <Box sx={{ textAlign: 'center', mb: 3 }}>
                        <PersonIcon sx={{ fontSize: 48, color: '#74a9cf', mb: 1 }} />
                        <Typography variant="h4" component="h1" gutterBottom sx={{ color: '#333', fontWeight: 600 }}>
                            Westmap
                        </Typography>
                        <Typography variant="body2" sx={{ color: '#666', mb: 1 }}>
                            Weather and Energy System Tracking and Modeling Analytics Platform
                        </Typography>
                        {/* <Typography variant="h6" sx={{ color: '#666', mb: 2 }}>
                            {isSignUp ? 'Create Account' : 'Welcome Back'}
                        </Typography> */}
                    </Box>

                    {error && (
                        <Alert severity="error" sx={{ mb: 2 }}>
                            {error}
                        </Alert>
                    )}

                    <Box component="form" onSubmit={handleSubmit}>
                        <TextField
                            fullWidth
                            label="Email Address"
                            type="email"
                            value={email}
                            onChange={(e) => setEmail(e.target.value)}
                            required
                            margin="normal"
                            variant="outlined"
                            sx={{
                                '& .MuiOutlinedInput-root': {
                                    '&:hover fieldset': {
                                        borderColor: '#74a9cf',
                                    },
                                    '&.Mui-focused fieldset': {
                                        borderColor: '#74a9cf',
                                    },
                                },
                                '& .MuiInputLabel-root.Mui-focused': {
                                    color: '#74a9cf',
                                },
                            }}
                        />

                        {isSignUp && (
                            <TextField
                                fullWidth
                                label="Display Name (Optional)"
                                type="text"
                                value={displayName}
                                onChange={(e) => setDisplayName(e.target.value)}
                                margin="normal"
                                variant="outlined"
                                sx={{
                                    '& .MuiOutlinedInput-root': {
                                        '&:hover fieldset': {
                                            borderColor: '#74a9cf',
                                        },
                                        '&.Mui-focused fieldset': {
                                            borderColor: '#74a9cf',
                                        },
                                    },
                                    '& .MuiInputLabel-root.Mui-focused': {
                                        color: '#74a9cf',
                                    },
                                }}
                            />
                        )}

                        <TextField
                            fullWidth
                            label="Password"
                            type="password"
                            value={password}
                            onChange={(e) => setPassword(e.target.value)}
                            required
                            margin="normal"
                            variant="outlined"
                            sx={{
                                '& .MuiOutlinedInput-root': {
                                    '&:hover fieldset': {
                                        borderColor: '#74a9cf',
                                    },
                                    '&.Mui-focused fieldset': {
                                        borderColor: '#74a9cf',
                                    },
                                },
                                '& .MuiInputLabel-root.Mui-focused': {
                                    color: '#74a9cf',
                                },
                            }}
                        />

                        {isSignUp && (
                            <TextField
                                fullWidth
                                label="Confirm Password"
                                type="password"
                                value={confirmPassword}
                                onChange={(e) => setConfirmPassword(e.target.value)}
                                required
                                margin="normal"
                                variant="outlined"
                                sx={{
                                    '& .MuiOutlinedInput-root': {
                                        '&:hover fieldset': {
                                            borderColor: '#74a9cf',
                                        },
                                        '&.Mui-focused fieldset': {
                                            borderColor: '#74a9cf',
                                        },
                                    },
                                    '& .MuiInputLabel-root.Mui-focused': {
                                        color: '#74a9cf',
                                    },
                                }}
                            />
                        )}

                        <Button
                            type="submit"
                            fullWidth
                            variant="contained"
                            disabled={loading}
                            startIcon={<LockIcon />}
                            sx={{
                                mt: 3,
                                mb: 2,
                                backgroundColor: '#74a9cf',
                                '&:hover': {
                                    backgroundColor: '#5a8bb5',
                                },
                                py: 1.5,
                                fontSize: '1.1rem',
                                fontWeight: 600
                            }}
                        >
                            {loading ? 'Please wait...' : (isSignUp ? 'Create Account' : 'Sign In')}
                        </Button>
                    </Box>

                    <Box sx={{ textAlign: 'center', mt: 2 }}>
                        <Typography variant="body2" sx={{ color: '#666' }}>
                            {isSignUp ? 'Already have an account?' : 'Don\'t have an account?'}
                            <Button
                                onClick={toggleMode}
                                sx={{
                                    ml: 1,
                                    color: '#74a9cf',
                                    textTransform: 'none',
                                    fontWeight: 600,
                                    '&:hover': {
                                        backgroundColor: 'rgba(116, 169, 207, 0.04)',
                                    }
                                }}
                            >
                                {isSignUp ? 'Sign In' : 'Sign Up'}
                            </Button>
                        </Typography>
                    </Box>
                </Paper>
            </Box>
        </ThemeProvider>
    );
};

export default Login;
