import React from 'react';
import { Box, Typography, Paper } from '@mui/material';

const ManishProject = () => {
    return (
        <Box sx={{
            width: '100%',
            height: '100vh',
            backgroundColor: '#f5f5f5',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center'
        }}>
            <Paper sx={{
                padding: 4,
                textAlign: 'center',
                maxWidth: 600
            }}>
                <Typography variant="h4" component="h1" gutterBottom>
                    Manish's Project (Placeholder)
                </Typography>
                <Typography variant="body1" color="textSecondary">
                    This placeholder should not be visible anymore.
                    The main map and chat widget UI is now served as Manish's project.
                </Typography>
            </Paper>
        </Box>
    );
};

export default ManishProject;
