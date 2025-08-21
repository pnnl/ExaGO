# Admin Setup Instructions

## Step 1: Complete Firebase Console Setup

### 1.1 Create Firestore Database
1. Go to your Firebase Console: https://console.firebase.google.com
2. Select your project
3. Navigate to **Firestore Database** in the left sidebar
4. Click **Create database**
5. Choose **Start in production mode**
6. Select your preferred location
7. Click **Done**

### 1.2 Set Up Firestore Security Rules
1. In Firestore Database, go to **Rules** tab
2. Replace the default rules with this code:

```javascript
rules_version = '2';
service cloud.firestore {
  match /databases/{database}/documents {
    // Users collection rules
    match /users/{userId} {
      // Users can read their own document
      allow read: if request.auth != null && request.auth.uid == userId;
      
      // Users can create their own document (only on registration)
      allow create: if request.auth != null && 
                   request.auth.uid == userId &&
                   request.resource.data.approved == false &&
                   request.resource.data.role == 'user';
      
      // Users can update their own profile (but not approval status or role)
      allow update: if request.auth != null && 
                   request.auth.uid == userId &&
                   request.resource.data.approved == resource.data.approved &&
                   request.resource.data.role == resource.data.role;
    }
    
    // Admin access to all user documents
    match /users/{userId} {
      allow read, write: if request.auth != null && 
                        exists(/databases/$(database)/documents/users/$(request.auth.uid)) &&
                        get(/databases/$(database)/documents/users/$(request.auth.uid)).data.role == 'admin' &&
                        get(/databases/$(database)/documents/users/$(request.auth.uid)).data.approved == true;
    }
  }
}
```

3. Click **Publish**

## Step 2: Create First Admin User

### Method 1: Register and Manually Promote (Recommended)

1. **Start your application**: `npm start`
2. **Register a new account** through the signup form
3. **Note down your Firebase Auth UID**:
   - Go to Firebase Console > Authentication > Users
   - Find your email and copy the User UID

4. **Create admin user document in Firestore**:
   - Go to Firebase Console > Firestore Database > Data tab
   - Click **Start collection**
   - Collection ID: `users`
   - Document ID: Paste your Firebase Auth UID
   - Add these fields:
     ```
     email: "your-email@example.com" (string)
     displayName: "Admin User" (string)
     approved: true (boolean)
     role: "admin" (string)
     createdAt: [current timestamp]
     emailVerified: true (boolean)
     uid: "your-firebase-auth-uid" (string)
     ```
   - Click **Save**

5. **Test admin access**:
   - Sign out and sign back in
   - Go to: `http://localhost:3000?admin=true`
   - You should see the admin dashboard

### Method 2: Direct Database Creation

If you prefer to create the admin user directly in Firestore:

1. Go to Firebase Console > Firestore Database > Data tab
2. Click **Start collection**
3. Collection ID: `users`
4. Document ID: Create any unique ID (or use a Firebase Auth UID)
5. Add the fields as shown above
6. Create the corresponding Firebase Auth user if needed

## Step 3: Test the System

### 3.1 Test Regular User Flow
1. **Register a new user**: Should create account but show "Pending Approval"
2. **Try to access app**: Should see pending approval screen
3. **Admin approves user**: Go to admin panel and approve the user
4. **User accesses app**: Should now be able to access the full application

### 3.2 Test Admin Functions
1. **Access admin panel**: `http://localhost:3000?admin=true`
2. **View user list**: See all registered users
3. **Approve users**: Test approval functionality
4. **Promote users**: Test promoting users to admin

## Step 4: URLs and Access

### For Users
- **Main App**: `http://localhost:3000`
- **Registration/Login**: Automatic redirect if not authenticated

### For Admins
- **Admin Dashboard**: `http://localhost:3000?admin=true`
- **Main App**: `http://localhost:3000` (admins can access both)

## Step 5: Deployment Considerations

When deploying to production:

1. **Update Firebase Rules**: Ensure production rules are properly set
2. **Environment Variables**: Set up production Firebase config
3. **Admin URLs**: Consider using a proper routing system for better admin URLs
4. **Security**: Consider additional security measures for admin access

## Troubleshooting

### Common Issues

1. **"Access Denied" on Admin Panel**
   - Check that your user document has `role: "admin"` and `approved: true`
   - Verify you're accessing the URL with `?admin=true`

2. **Users Stuck in Pending State**
   - Check Firebase Console for user documents
   - Verify Firestore rules are properly set

3. **Firestore Permission Errors**
   - Double-check Firestore security rules
   - Ensure rules are published

4. **Admin Dashboard Not Loading**
   - Check browser console for errors
   - Verify all Firebase imports are correct

## Next Steps

After setup is complete:
1. Train admins on how to use the approval system
2. Consider adding email notifications for user approval
3. Set up monitoring for user registration patterns
4. Consider implementing bulk approval features if needed

---

**Important**: Keep your Firebase configuration and admin credentials secure. Never share admin access credentials publicly.
