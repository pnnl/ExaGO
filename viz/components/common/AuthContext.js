import React, { createContext, useContext, useState, useEffect } from 'react';
import {
    signInWithEmailAndPassword,
    createUserWithEmailAndPassword,
    signOut,
    onAuthStateChanged,
    updateProfile
} from 'firebase/auth';
import {
    doc,
    setDoc,
    getDoc,
    updateDoc,
    collection,
    query,
    where,
    onSnapshot
} from 'firebase/firestore';
import { auth, db } from '../../firebase';

const AuthContext = createContext();

export function useAuth() {
    return useContext(AuthContext);
}

export function AuthProvider({ children }) {
    const [currentUser, setCurrentUser] = useState(null);
    const [userProfile, setUserProfile] = useState(null);
    const [loading, setLoading] = useState(true);

    function login(email, password) {
        return signInWithEmailAndPassword(auth, email, password);
    }

    async function signup(email, password, displayName = '') {
        const { user } = await createUserWithEmailAndPassword(auth, email, password);

        // Update profile with display name if provided
        if (displayName) {
            await updateProfile(user, { displayName });
        }

        // Create user document in Firestore
        await setDoc(doc(db, 'users', user.uid), {
            email: user.email,
            displayName: displayName || '',
            approved: false,
            role: 'user',
            createdAt: new Date(),
            emailVerified: user.emailVerified,
            uid: user.uid
        });

        return user;
    }

    function logout() {
        return signOut(auth);
    }

    // Admin functions
    async function approveUser(userId) {
        await updateDoc(doc(db, 'users', userId), {
            approved: true
        });
    }

    async function rejectUser(userId) {
        await updateDoc(doc(db, 'users', userId), {
            approved: false
        });
    }

    async function promoteToAdmin(userId) {
        await updateDoc(doc(db, 'users', userId), {
            role: 'admin'
        });
    }

    async function getUserProfile(userId) {
        const userDoc = await getDoc(doc(db, 'users', userId));
        return userDoc.exists() ? userDoc.data() : null;
    }

    useEffect(() => {
        const unsubscribe = onAuthStateChanged(auth, async (user) => {
            setCurrentUser(user);

            if (user) {
                // Get user profile from Firestore
                try {
                    const profile = await getUserProfile(user.uid);
                    setUserProfile(profile);
                } catch (error) {
                    console.error('Error fetching user profile:', error);
                    setUserProfile(null);
                }
            } else {
                setUserProfile(null);
            }

            setLoading(false);
        });

        return unsubscribe;
    }, []);

    const value = {
        currentUser,
        userProfile,
        login,
        signup,
        logout,
        approveUser,
        rejectUser,
        promoteToAdmin,
        getUserProfile
    };

    return (
        <AuthContext.Provider value={value}>
            {!loading && children}
        </AuthContext.Provider>
    );
}
