# Westmap - Authentication System

## Overview
This project implements Firebase Authentication for the Westmap application with a clean, Material-UI based interface that matches the app's existing design theme.

## Features

### ✅ Implemented Features
- **Email/Password Authentication**: Users can sign up and log in with email and password
- **Protected Routes**: Main visualization is only accessible after authentication
- **Persistent Sessions**: Users remain logged in across browser sessions
- **Responsive Design**: Auth components work on desktop and mobile devices
- **Custom Styling**: Matches the app's color scheme (`#74a9cf` primary color)

### 🔧 Technical Implementation
- **Firebase SDK v12**: Latest Firebase authentication features
- **Material-UI Components**: Consistent design system
- **React Context**: Centralized auth state management
- **Custom Theme**: Branded colors and styling
- **Proper Folder Structure**: Organized auth components

## File Structure

```
viz/
├── components/
│   └── auth/
│       ├── index.js          # Barrel exports
│       ├── AuthContext.js    # Authentication context & hooks
│       ├── Login.js          # Login/Signup form component
│       ├── ProtectedRoute.js # Route protection wrapper
│       ├── Header.js         # User info & logout header
│       └── theme.js          # Custom Material-UI theme
├── firebase.js               # Firebase configuration
├── .env                      # Environment variables
└── app.js                    # Main app with auth integration
```

## Environment Variables

The following environment variables are required in `.env`:

```env
# Firebase Configuration
FIREBASE_API_KEY=your_api_key
FIREBASE_AUTH_DOMAIN=your_project.firebaseapp.com
FIREBASE_PROJECT_ID=your_project_id
FIREBASE_STORAGE_BUCKET=your_project.firebasestorage.app
FIREBASE_MESSAGING_SENDER_ID=your_sender_id
FIREBASE_APP_ID=your_app_id
FIREBASE_MEASUREMENT_ID=your_measurement_id
```

## Usage

### For Users
1. **Access the App**: Navigate to the application URL
2. **Sign Up**: Click "Sign Up" to create a new account with email/password
3. **Sign In**: Use existing credentials with email/password
4. **Access Visualization**: Authenticated users can access the full visualization
5. **Logout**: Click the logout button in the top-right header

### For Developers

#### Authentication Context
```javascript
import { useAuth } from './components/auth';

function YourComponent() {
  const { currentUser, login, signup, logout } = useAuth();
  // Use auth methods and current user state
}
```

#### Protected Routes
```javascript
import { ProtectedRoute } from './components/auth';

function App() {
  return (
    <ProtectedRoute>
      <YourProtectedComponent />
    </ProtectedRoute>
  );
}
```

## Security Features

- **Environment Variables**: Sensitive Firebase config stored in `.env`
- **Client-side Validation**: Form validation before submission
- **Firebase Security Rules**: Server-side authentication enforcement
- **Automatic Token Refresh**: Firebase handles token lifecycle
- **Error Handling**: User-friendly error messages

## Styling

The authentication components use a custom Material-UI theme that matches the main application:

- **Primary Color**: `#74a9cf` (matches chat widget and other UI elements)
- **Secondary Color**: `#dc3545` (for logout and error states)  
- **Background**: Gradient from `#74a9cf` to `#a8c8e1`
- **Modern Design**: Cards, shadows, and smooth transitions

## Testing

To test the authentication system:

1. **Start Development Server**: `npm run start`
2. **Test Signup**: Create a new account with email/password
3. **Test Login**: Sign in with existing credentials
4. **Test Logout**: Verify logout functionality
5. **Test Protection**: Verify unauthenticated users see login form

## Troubleshooting

### Common Issues

1. **Firebase Config Error**
   - Check `.env` file has all required variables
   - Verify Firebase project is properly configured

2. **Build Errors**
   - Ensure all dependencies are installed: `npm install`
   - Check for missing Material-UI components

### Error Messages

The system provides user-friendly error messages for:
- Invalid email format
- Password too weak
- Account doesn't exist
- Network connection issues
- Firebase service errors

## Future Enhancements

Potential improvements for the auth system:

- **Email Verification**: Require email verification for new accounts
- **Password Reset**: "Forgot Password" functionality
- **Profile Management**: User profile editing
- **Role-based Access**: Different user permissions
- **Social Logins**: Facebook, GitHub, etc.
- **Multi-factor Authentication**: Enhanced security
- **User Analytics**: Login tracking and metrics

## Dependencies

Key dependencies for the auth system:

```json
{
  "firebase": "^12.0.0",
  "@mui/material": "^5.6.4",
  "@mui/icons-material": "^5.6.2",
  "@emotion/react": "^11.9.0",
  "@emotion/styled": "^11.8.1"
}
```

## Support

For issues or questions about the authentication system:

1. Check the console for error messages
2. Verify Firebase configuration
3. Review this documentation
4. Check Firebase Console for user management

---

**Note**: This authentication system is designed to integrate seamlessly with the existing Westmap application while maintaining security best practices and a consistent user experience.
