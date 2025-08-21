import React from 'react';
import { useAuth } from './AuthContext';
import { Box, Button, Typography, Avatar, Paper, ThemeProvider, Chip } from '@mui/material';
import LogoutIcon from '@mui/icons-material/Logout';
import PersonIcon from '@mui/icons-material/Person';
import AdminPanelSettingsIcon from '@mui/icons-material/AdminPanelSettings';
import { authTheme } from './theme';

function Header() {
    const { currentUser, userProfile, logout } = useAuth();

    async function handleLogout() {
        try {
            await logout();
        } catch (error) {
            console.error('Failed to log out:', error);
        }
    }

    if (!currentUser) return null;

    return (
        <ThemeProvider theme={authTheme}>
            <Paper
                elevation={3}
                sx={{
                    position: 'fixed',
                    top: 0,
                    left: '50vw',
                    transform: 'translateX(-50%)',
                    zIndex: 1200,
                    backgroundColor: 'rgba(255, 255, 255, 0.95)',
                    backdropFilter: 'blur(10px)',
                    padding: '12px 20px',
                    borderBottomLeftRadius: 12,
                    display: 'flex',
                    alignItems: 'center',
                    gap: 2,
                    minWidth: 250
                }}
            >
                <Avatar
                    sx={{
                        width: 32,
                        height: 32,
                        backgroundColor: '#74a9cf',
                        fontSize: '14px'
                    }}
                    src={currentUser.photoURL}
                >
                    {!currentUser.photoURL && <PersonIcon sx={{ fontSize: 18 }} />}
                </Avatar>

                <Box sx={{ flex: 1 }}>
                    <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, mb: 0.5 }}>
                        <Typography
                            variant="body2"
                            sx={{
                                fontWeight: 600,
                                color: '#333',
                                fontSize: '14px'
                            }}
                        >
                            {currentUser.displayName || userProfile?.displayName || 'User'}
                        </Typography>
                        {userProfile?.role === 'admin' && (
                            <Chip
                                icon={<AdminPanelSettingsIcon sx={{ fontSize: 12 }} />}
                                label="Admin"
                                size="small"
                                color="secondary"
                                sx={{ height: 16, fontSize: '10px' }}
                            />
                        )}
                    </Box>
                    <Typography
                        variant="caption"
                        sx={{
                            color: '#666',
                            fontSize: '12px',
                            display: 'block',
                            lineHeight: 1
                        }}
                    >
                        {currentUser.email}
                    </Typography>
                </Box>

                {/* Project Navigation */}
                <Box sx={{ display: 'flex', gap: 1, alignItems: 'center' }}>
                    {/* Manish's Project */}
                    <Button
                        onClick={() => {
                            window.history.pushState({}, '', '/manish');
                            window.location.reload();
                        }}
                        size="small"
                        sx={{
                            color: '#293842ff',
                            fontSize: '10px',
                            textTransform: 'none',
                            minWidth: 'auto',
                            padding: '2px 6px',
                            '&:hover': {
                                backgroundColor: 'rgba(41, 56, 66, 1)',
                                color: 'white',
                            }
                        }}
                    >
                        Manish's Project
                    </Button>

                    {/* Amin's Project */}
                    <Button
                        onClick={() => {
                            window.history.pushState({}, '', '/amin');
                            window.location.reload();
                        }}
                        size="small"
                        sx={{
                            color: '#293842ff',
                            fontSize: '10px',
                            textTransform: 'none',
                            minWidth: 'auto',
                            padding: '2px 6px',
                            '&:hover': {
                                backgroundColor: 'rgba(41, 56, 66, 1)',
                                color: 'white',
                            }
                        }}
                    >
                        Amin's Project
                    </Button>

                    {/* Admin Panel - only for admin users */}
                    {userProfile?.role === 'admin' && (
                        <Button
                            onClick={() => {
                                window.history.pushState({}, '', '/admin');
                                window.location.reload();
                            }}
                            size="small"
                            sx={{
                                color: '#9c27b0',
                                fontSize: '10px',
                                textTransform: 'none',
                                minWidth: 'auto',
                                padding: '2px 6px',
                                border: '1px solid #9c27b0',
                                '&:hover': {
                                    backgroundColor: '#9c27b0',
                                    color: 'white',
                                }
                            }}
                        >
                            Admin Panel
                        </Button>
                    )}
                </Box>

                <Button
                    onClick={handleLogout}
                    size="small"
                    startIcon={<LogoutIcon sx={{ fontSize: 16 }} />}
                    sx={{
                        padding: '6px 12px',
                        backgroundColor: '#dc3545',
                        color: 'white',
                        fontSize: '12px',
                        fontWeight: 500,
                        borderRadius: 2,
                        minWidth: 'auto',
                        '&:hover': {
                            backgroundColor: '#c82333',
                        }
                    }}
                >
                    Logout
                </Button>
            </Paper>
        </ThemeProvider>
    );
}

export default Header;
