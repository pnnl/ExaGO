import React from 'react';
import { useAuth } from './AuthContext';
import Login from './Login';
import PendingApproval from './PendingApproval';

function ProtectedRoute({ children, adminOnly = false }) {
    const { currentUser, userProfile } = useAuth();

    // Not logged in - show login
    if (!currentUser) {
        return <Login />;
    }

    // Logged in but no profile data yet - show loading
    if (!userProfile) {
        return <div>Loading...</div>;
    }

    // User not approved - show pending approval screen
    if (!userProfile.approved) {
        return <PendingApproval />;
    }

    // Admin-only route but user is not admin
    if (adminOnly && userProfile.role !== 'admin') {
        return (
            <div style={{ padding: '20px', textAlign: 'center' }}>
                <h2>Access Denied</h2>
                <p>You need admin privileges to access this page.</p>
            </div>
        );
    }

    // User is approved - show protected content
    return children;
}

export default ProtectedRoute;
