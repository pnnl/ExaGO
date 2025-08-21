import React, { useState, useEffect } from 'react';
import {
    Box,
    Paper,
    Typography,
    Button,
    Table,
    TableBody,
    TableCell,
    TableContainer,
    TableHead,
    TableRow,
    Chip,
    IconButton,
    Dialog,
    DialogTitle,
    DialogContent,
    DialogActions,
    Alert,
    Tabs,
    Tab,
    ThemeProvider,
    Tooltip,
    Card,
    CardContent,
    Grid
} from '@mui/material';
import {
    Check as CheckIcon,
    Close as CloseIcon,
    AdminPanelSettings as AdminIcon,
    People as PeopleIcon,
    HourglassEmpty as PendingIcon,
    CheckCircle as ApprovedIcon,
    Cancel as RejectedIcon,
    Delete as DeleteIcon,
    Sort as SortIcon,
    Dashboard as DashboardIcon,
    Work as WorkIcon
} from '@mui/icons-material';
import { collection, onSnapshot, query, orderBy, where, deleteDoc, doc } from 'firebase/firestore';
import { useAuth } from './AuthContext';
import { db } from '../../firebase';
import { authTheme } from './theme';

const AdminDashboard = () => {
    const { currentUser, userProfile, approveUser, rejectUser, promoteToAdmin } = useAuth();
    const [users, setUsers] = useState([]);
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState('');
    const [selectedUser, setSelectedUser] = useState(null);
    const [actionDialog, setActionDialog] = useState({ open: false, action: '', user: null });
    const [tabValue, setTabValue] = useState(0);
    const [sortBy, setSortBy] = useState('createdAt');
    const [sortOrder, setSortOrder] = useState('desc');
    const [filterStatus, setFilterStatus] = useState('all');

    // Statistics
    const totalUsers = users.length;
    const pendingUsers = users.filter(user => !user.approved).length;
    const approvedUsers = users.filter(user => user.approved && user.role !== 'admin').length;
    const adminUsers = users.filter(user => user.role === 'admin').length;
    const rejectedUsers = users.filter(user => user.approved === false && user.rejected === true).length;

    useEffect(() => {
        if (!currentUser || !userProfile?.approved || userProfile?.role !== 'admin') {
            return;
        }

        const usersQuery = query(
            collection(db, 'users'),
            orderBy('createdAt', 'desc')
        );

        const unsubscribe = onSnapshot(
            usersQuery,
            (snapshot) => {
                const usersData = snapshot.docs.map(doc => ({
                    id: doc.id,
                    ...doc.data()
                }));
                setUsers(usersData);
                setLoading(false);
            },
            (error) => {
                console.error('Error fetching users:', error);
                setError('Failed to load users');
                setLoading(false);
            }
        );

        return unsubscribe;
    }, [currentUser, userProfile]);

    // Admin functions
    async function deleteUser(userId) {
        await deleteDoc(doc(db, 'users', userId));
    }

    const handleAction = async (action, user) => {
        try {
            setError('');

            switch (action) {
                case 'approve':
                    await approveUser(user.id);
                    break;
                case 'reject':
                    await rejectUser(user.id);
                    break;
                case 'promote':
                    await promoteToAdmin(user.id);
                    break;
                case 'delete':
                    await deleteUser(user.id);
                    break;
                default:
                    throw new Error('Invalid action');
            }

            setActionDialog({ open: false, action: '', user: null });
        } catch (error) {
            console.error('Error performing action:', error);
            setError(`Failed to ${action} user: ${error.message}`);
        }
    };

    const getFilteredUsers = () => {
        let filtered = users;

        // Apply tab filter
        switch (tabValue) {
            case 0: // All
                filtered = users;
                break;
            case 1: // Pending
                filtered = users.filter(user => !user.approved);
                break;
            case 2: // Approved
                filtered = users.filter(user => user.approved && user.role !== 'admin');
                break;
            case 3: // Admins
                filtered = users.filter(user => user.role === 'admin');
                break;
            default:
                filtered = users;
        }

        // Apply sorting
        filtered.sort((a, b) => {
            let aValue, bValue;

            switch (sortBy) {
                case 'email':
                    aValue = a.email || '';
                    bValue = b.email || '';
                    break;
                case 'displayName':
                    aValue = a.displayName || '';
                    bValue = b.displayName || '';
                    break;
                case 'role':
                    aValue = a.role || '';
                    bValue = b.role || '';
                    break;
                case 'approved':
                    aValue = a.approved ? 1 : 0;
                    bValue = b.approved ? 1 : 0;
                    break;
                case 'createdAt':
                default:
                    aValue = a.createdAt?.toDate?.() || new Date(0);
                    bValue = b.createdAt?.toDate?.() || new Date(0);
                    break;
            }

            if (sortOrder === 'asc') {
                return aValue > bValue ? 1 : -1;
            } else {
                return aValue < bValue ? 1 : -1;
            }
        });

        return filtered;
    };

    const handleSort = (column) => {
        if (sortBy === column) {
            setSortOrder(sortOrder === 'asc' ? 'desc' : 'asc');
        } else {
            setSortBy(column);
            setSortOrder('asc');
        }
    };

    const getUserStatus = (user) => {
        if (user.role === 'admin') {
            return { label: 'Admin', color: 'secondary' };
        }
        if (user.approved) {
            return { label: 'Approved', color: 'success' };
        }
        return { label: 'Pending', color: 'warning' };
    };

    if (!userProfile?.approved || userProfile?.role !== 'admin') {
        return (
            <Box sx={{ p: 3, textAlign: 'center' }}>
                <Typography variant="h6" color="error">
                    Access Denied: Admin privileges required
                </Typography>
            </Box>
        );
    }

    return (
        <ThemeProvider theme={authTheme}>
            <Box sx={{ p: 3, backgroundColor: '#f5f5f5', minHeight: '100vh' }}>
                <Typography variant="h4" gutterBottom sx={{ color: '#333', fontWeight: 600 }}>
                    Admin Dashboard
                </Typography>

                {error && (
                    <Alert severity="error" sx={{ mb: 3 }}>
                        {error}
                    </Alert>
                )}

                {/* Statistics Cards */}
                <Grid container spacing={3} sx={{ mb: 4 }}>
                    <Grid item xs={12} sm={6} md={3}>
                        <Card>
                            <CardContent sx={{ textAlign: 'center' }}>
                                <PeopleIcon sx={{ fontSize: 40, color: '#74a9cf', mb: 1 }} />
                                <Typography variant="h4" sx={{ fontWeight: 600 }}>
                                    {totalUsers}
                                </Typography>
                                <Typography variant="body2" color="textSecondary">
                                    Total Users
                                </Typography>
                            </CardContent>
                        </Card>
                    </Grid>
                    <Grid item xs={12} sm={6} md={3}>
                        <Card>
                            <CardContent sx={{ textAlign: 'center' }}>
                                <PendingIcon sx={{ fontSize: 40, color: '#ff9800', mb: 1 }} />
                                <Typography variant="h4" sx={{ fontWeight: 600 }}>
                                    {pendingUsers}
                                </Typography>
                                <Typography variant="body2" color="textSecondary">
                                    Pending Approval
                                </Typography>
                            </CardContent>
                        </Card>
                    </Grid>
                    <Grid item xs={12} sm={6} md={3}>
                        <Card>
                            <CardContent sx={{ textAlign: 'center' }}>
                                <ApprovedIcon sx={{ fontSize: 40, color: '#4caf50', mb: 1 }} />
                                <Typography variant="h4" sx={{ fontWeight: 600 }}>
                                    {approvedUsers}
                                </Typography>
                                <Typography variant="body2" color="textSecondary">
                                    Approved Users
                                </Typography>
                            </CardContent>
                        </Card>
                    </Grid>
                    <Grid item xs={12} sm={6} md={3}>
                        <Card>
                            <CardContent sx={{ textAlign: 'center' }}>
                                <AdminIcon sx={{ fontSize: 40, color: '#9c27b0', mb: 1 }} />
                                <Typography variant="h4" sx={{ fontWeight: 600 }}>
                                    {adminUsers}
                                </Typography>
                                <Typography variant="body2" color="textSecondary">
                                    Administrators
                                </Typography>
                            </CardContent>
                        </Card>
                    </Grid>
                </Grid>

                {/* User Filter Tabs */}
                <Paper sx={{ mb: 3 }}>
                    <Box sx={{ p: 2, borderBottom: 1, borderColor: 'divider' }}>
                        <Typography variant="h6" sx={{ mb: 2 }}>User Management</Typography>
                        <Box sx={{ display: 'flex', gap: 2, alignItems: 'center', flexWrap: 'wrap' }}>
                            <Typography variant="body2" sx={{ fontWeight: 600 }}>
                                Sort by:
                            </Typography>
                            <Button
                                size="small"
                                variant={sortBy === 'createdAt' ? 'contained' : 'outlined'}
                                onClick={() => handleSort('createdAt')}
                                sx={{ textTransform: 'none' }}
                            >
                                Date {sortBy === 'createdAt' && (sortOrder === 'asc' ? '↑' : '↓')}
                            </Button>
                            <Button
                                size="small"
                                variant={sortBy === 'email' ? 'contained' : 'outlined'}
                                onClick={() => handleSort('email')}
                                sx={{ textTransform: 'none' }}
                            >
                                Email {sortBy === 'email' && (sortOrder === 'asc' ? '↑' : '↓')}
                            </Button>
                            <Button
                                size="small"
                                variant={sortBy === 'role' ? 'contained' : 'outlined'}
                                onClick={() => handleSort('role')}
                                sx={{ textTransform: 'none' }}
                            >
                                Role {sortBy === 'role' && (sortOrder === 'asc' ? '↑' : '↓')}
                            </Button>
                            <Button
                                size="small"
                                variant={sortBy === 'approved' ? 'contained' : 'outlined'}
                                onClick={() => handleSort('approved')}
                                sx={{ textTransform: 'none' }}
                            >
                                Status {sortBy === 'approved' && (sortOrder === 'asc' ? '↑' : '↓')}
                            </Button>
                        </Box>
                    </Box>
                    <Tabs
                        value={tabValue}
                        onChange={(e, newValue) => setTabValue(newValue)}
                        sx={{ borderBottom: 1, borderColor: 'divider' }}
                    >
                        <Tab label={`All Users (${totalUsers})`} />
                        <Tab label={`Pending Users (${pendingUsers})`} />
                        <Tab label={`Approved Users (${approvedUsers})`} />
                        <Tab label={`Admins (${adminUsers})`} />
                    </Tabs>
                </Paper>

                {/* Users Table */}
                <TableContainer component={Paper}>
                    <Table>
                        <TableHead>
                            <TableRow>
                                <TableCell>
                                    <Box
                                        sx={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}
                                        onClick={() => handleSort('email')}
                                    >
                                        <strong>Email</strong>
                                        <SortIcon sx={{ ml: 1, fontSize: 16 }} />
                                    </Box>
                                </TableCell>
                                <TableCell>
                                    <Box
                                        sx={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}
                                        onClick={() => handleSort('displayName')}
                                    >
                                        <strong>Display Name</strong>
                                        <SortIcon sx={{ ml: 1, fontSize: 16 }} />
                                    </Box>
                                </TableCell>
                                <TableCell>
                                    <Box
                                        sx={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}
                                        onClick={() => handleSort('approved')}
                                    >
                                        <strong>Status</strong>
                                        <SortIcon sx={{ ml: 1, fontSize: 16 }} />
                                    </Box>
                                </TableCell>
                                <TableCell>
                                    <Box
                                        sx={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}
                                        onClick={() => handleSort('role')}
                                    >
                                        <strong>Role</strong>
                                        <SortIcon sx={{ ml: 1, fontSize: 16 }} />
                                    </Box>
                                </TableCell>
                                <TableCell>
                                    <Box
                                        sx={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}
                                        onClick={() => handleSort('createdAt')}
                                    >
                                        <strong>Created At</strong>
                                        <SortIcon sx={{ ml: 1, fontSize: 16 }} />
                                    </Box>
                                </TableCell>
                                <TableCell><strong>Actions</strong></TableCell>
                            </TableRow>
                        </TableHead>
                        <TableBody>
                            {getFilteredUsers().map((user) => {
                                const status = getUserStatus(user);
                                return (
                                    <TableRow key={user.id}>
                                        <TableCell>{user.email}</TableCell>
                                        <TableCell>{user.displayName || '-'}</TableCell>
                                        <TableCell>
                                            <Chip
                                                label={status.label}
                                                color={status.color}
                                                size="small"
                                            />
                                        </TableCell>
                                        <TableCell>
                                            <Chip
                                                label={user.role}
                                                variant="outlined"
                                                size="small"
                                            />
                                        </TableCell>
                                        <TableCell>
                                            {user.createdAt?.toDate?.()?.toLocaleDateString() || 'N/A'}
                                        </TableCell>
                                        <TableCell>
                                            <Box sx={{ display: 'flex', gap: 1, flexWrap: 'wrap' }}>
                                                {/* Approve button for pending users */}
                                                {!user.approved && (
                                                    <Tooltip title="Approve User">
                                                        <IconButton
                                                            color="success"
                                                            size="small"
                                                            onClick={() => setActionDialog({
                                                                open: true,
                                                                action: 'approve',
                                                                user
                                                            })}
                                                        >
                                                            <CheckIcon />
                                                        </IconButton>
                                                    </Tooltip>
                                                )}

                                                {/* Reject button for pending or approved users (not admins) */}
                                                {user.role !== 'admin' && user.id !== currentUser.uid && (
                                                    <Tooltip title="Reject User">
                                                        <IconButton
                                                            color="error"
                                                            size="small"
                                                            onClick={() => setActionDialog({
                                                                open: true,
                                                                action: 'reject',
                                                                user
                                                            })}
                                                        >
                                                            <CloseIcon />
                                                        </IconButton>
                                                    </Tooltip>
                                                )}

                                                {/* Promote to Admin button for approved users */}
                                                {user.approved && user.role !== 'admin' && user.id !== currentUser.uid && (
                                                    <Tooltip title="Promote to Admin">
                                                        <IconButton
                                                            color="secondary"
                                                            size="small"
                                                            onClick={() => setActionDialog({
                                                                open: true,
                                                                action: 'promote',
                                                                user
                                                            })}
                                                        >
                                                            <AdminIcon />
                                                        </IconButton>
                                                    </Tooltip>
                                                )}

                                                {/* Delete button for all users except current user */}
                                                {user.id !== currentUser.uid && (
                                                    <Tooltip title="Delete User">
                                                        <IconButton
                                                            color="error"
                                                            size="small"
                                                            sx={{ opacity: 0.7 }}
                                                            onClick={() => setActionDialog({
                                                                open: true,
                                                                action: 'delete',
                                                                user
                                                            })}
                                                        >
                                                            <DeleteIcon />
                                                        </IconButton>
                                                    </Tooltip>
                                                )}

                                                {/* Current user indicator */}
                                                {user.id === currentUser.uid && (
                                                    <Chip
                                                        label="You"
                                                        color="primary"
                                                        size="small"
                                                        variant="outlined"
                                                    />
                                                )}
                                            </Box>
                                        </TableCell>
                                    </TableRow>
                                );
                            })}
                            {getFilteredUsers().length === 0 && (
                                <TableRow>
                                    <TableCell colSpan={6} sx={{ textAlign: 'center', py: 4 }}>
                                        <Typography variant="body2" color="textSecondary">
                                            No users found
                                        </Typography>
                                    </TableCell>
                                </TableRow>
                            )}
                        </TableBody>
                    </Table>
                </TableContainer>

                {/* Action Confirmation Dialog */}
                <Dialog
                    open={actionDialog.open}
                    onClose={() => setActionDialog({ open: false, action: '', user: null })}
                >
                    <DialogTitle>
                        Confirm {actionDialog.action?.charAt(0).toUpperCase() + actionDialog.action?.slice(1)}
                    </DialogTitle>
                    <DialogContent>
                        <Typography>
                            Are you sure you want to {actionDialog.action} user "{actionDialog.user?.email}"?
                            {actionDialog.action === 'delete' && (
                                <Box sx={{ mt: 2, p: 2, backgroundColor: '#ffebee', borderRadius: 1 }}>
                                    <Typography variant="body2" color="error" sx={{ fontWeight: 600 }}>
                                        ⚠️ Warning: This action cannot be undone. The user will be permanently deleted from both the database and Firebase Authentication.
                                    </Typography>
                                </Box>
                            )}
                        </Typography>
                    </DialogContent>
                    <DialogActions>
                        <Button
                            onClick={() => setActionDialog({ open: false, action: '', user: null })}
                        >
                            Cancel
                        </Button>
                        <Button
                            onClick={() => handleAction(actionDialog.action, actionDialog.user)}
                            color={actionDialog.action === 'delete' ? 'error' : 'primary'}
                            variant="contained"
                        >
                            {actionDialog.action === 'delete' ? 'Delete Permanently' : 'Confirm'}
                        </Button>
                    </DialogActions>
                </Dialog>
            </Box>
        </ThemeProvider>
    );
};

export default AdminDashboard;
