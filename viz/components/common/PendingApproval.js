import React from 'react';
import { Box, Paper, Typography, Button, ThemeProvider } from '@mui/material';
import HourglassEmptyIcon from '@mui/icons-material/HourglassEmpty';
import LogoutIcon from '@mui/icons-material/Logout';
import { useAuth } from './AuthContext';
import { authTheme } from './theme';

const PendingApproval = () => {
    const { currentUser, logout } = useAuth();

    const handleLogout = async () => {
        try {
            await logout();
        } catch (error) {
            console.error('Failed to log out:', error);
        }
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
                        maxWidth: 500,
                        backgroundColor: 'rgba(255, 255, 255, 0.95)',
                        backdropFilter: 'blur(10px)',
                        textAlign: 'center'
                    }}
                >
                    <HourglassEmptyIcon
                        sx={{
                            fontSize: 64,
                            color: '#74a9cf',
                            mb: 2
                        }}
                    />

                    <Typography
                        variant="h4"
                        component="h1"
                        gutterBottom
                        sx={{
                            color: '#333',
                            fontWeight: 600,
                            mb: 1
                        }}
                    >
                        Account Pending Approval
                    </Typography>

                    <Typography
                        variant="h6"
                        sx={{
                            color: '#666',
                            mb: 3,
                            fontWeight: 400
                        }}
                    >
                        Welcome to Westmap!
                    </Typography>

                    <Box
                        sx={{
                            backgroundColor: 'rgba(116, 169, 207, 0.1)',
                            padding: 3,
                            borderRadius: 2,
                            mb: 3
                        }}
                    >
                        <Typography
                            variant="body1"
                            sx={{
                                color: '#555',
                                lineHeight: 1.6,
                                mb: 2
                            }}
                        >
                            Your account has been successfully created! However, it requires admin approval before you can access the Westmap platform.
                        </Typography>

                        <Typography
                            variant="body2"
                            sx={{
                                color: '#666',
                                lineHeight: 1.5
                            }}
                        >
                            <strong>Account Details:</strong><br />
                            Email: {currentUser?.email}<br />
                            Status: Pending Admin Approval
                        </Typography>
                    </Box>

                    <Typography
                        variant="body2"
                        sx={{
                            color: '#777',
                            mb: 3,
                            fontStyle: 'italic'
                        }}
                    >
                        You will receive email notification once your account is approved. Please check your email regularly.
                    </Typography>

                    <Button
                        onClick={handleLogout}
                        variant="outlined"
                        startIcon={<LogoutIcon />}
                        sx={{
                            borderColor: '#74a9cf',
                            color: '#74a9cf',
                            '&:hover': {
                                backgroundColor: 'rgba(116, 169, 207, 0.04)',
                                borderColor: '#5a8bb5',
                            },
                            py: 1,
                            px: 3,
                            fontWeight: 600
                        }}
                    >
                        Sign Out
                    </Button>
                </Paper>
            </Box>
        </ThemeProvider>
    );
};

export default PendingApproval;
