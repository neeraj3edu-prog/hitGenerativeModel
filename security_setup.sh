#!/bin/bash

################################################################################
# Security Setup Script for OAM IITK Streamlit Application
# 
# This script sets up security measures for the Streamlit application including:
# - Nginx reverse proxy with SSL
# - Firewall configuration
# - Password hashing utilities
# - Session management
# 
# Usage: 
#   chmod +x security_setup.sh
#   sudo ./security_setup.sh
################################################################################

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
APP_DIR="$HOME/oam_iitk"
STREAMLIT_PORT=8501
DOMAIN=""  # Leave empty if using IP only

################################################################################
# Helper Functions
################################################################################

print_header() {
    echo -e "\n${GREEN}========================================${NC}"
    echo -e "${GREEN}$1${NC}"
    echo -e "${GREEN}========================================${NC}\n"
}

print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

check_root() {
    if [[ $EUID -ne 0 ]]; then
        print_error "This script must be run as root (use sudo)"
        exit 1
    fi
}

################################################################################
# Step 1: Install Required Packages
################################################################################

install_packages() {
    print_header "Step 1: Installing Required Packages"
    
    print_info "Updating package list..."
    apt update
    
    print_info "Installing Nginx, Certbot, and dependencies..."
    apt install -y nginx certbot python3-certbot-nginx ufw
    
    print_info "Packages installed successfully!"
}

################################################################################
# Step 2: Configure Firewall (UFW)
################################################################################

configure_firewall() {
    print_header "Step 2: Configuring Firewall"
    
    print_info "Setting up UFW firewall rules..."
    
    # Allow SSH
    ufw allow 22/tcp
    print_info "Allowed SSH (port 22)"
    
    # Allow HTTP and HTTPS
    ufw allow 80/tcp
    ufw allow 443/tcp
    print_info "Allowed HTTP (port 80) and HTTPS (port 443)"
    
    # Block direct access to Streamlit port from outside
    # (only allow from localhost)
    print_info "Blocking direct external access to Streamlit port $STREAMLIT_PORT"
    
    # Enable firewall
    print_warning "Enabling UFW firewall..."
    echo "y" | ufw enable
    
    ufw status
    print_info "Firewall configured successfully!"
}

################################################################################
# Step 3: Configure Nginx Reverse Proxy
################################################################################

configure_nginx() {
    print_header "Step 3: Configuring Nginx Reverse Proxy"
    
    # Get server IP or domain
    if [ -z "$DOMAIN" ]; then
        SERVER_IP=$(curl -s ifconfig.me)
        SERVER_NAME=$SERVER_IP
        print_info "Using IP address: $SERVER_IP"
    else
        SERVER_NAME=$DOMAIN
        print_info "Using domain: $DOMAIN"
    fi
    
    # Create Nginx configuration
    print_info "Creating Nginx configuration..."
    
    cat > /etc/nginx/sites-available/streamlit << EOF
server {
    listen 80;
    server_name $SERVER_NAME;

    # Increase buffer sizes for Streamlit
    client_max_body_size 100M;
    proxy_buffering off;

    # Security headers
    add_header X-Frame-Options "SAMEORIGIN" always;
    add_header X-Content-Type-Options "nosniff" always;
    add_header X-XSS-Protection "1; mode=block" always;
    add_header Referrer-Policy "no-referrer-when-downgrade" always;

    # Main application
    location / {
        proxy_pass http://localhost:$STREAMLIT_PORT;
        proxy_http_version 1.1;
        proxy_set_header Upgrade \$http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host \$host;
        proxy_set_header X-Real-IP \$remote_addr;
        proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto \$scheme;
        proxy_read_timeout 86400;
    }

    # WebSocket support for Streamlit
    location /_stcore/stream {
        proxy_pass http://localhost:$STREAMLIT_PORT/_stcore/stream;
        proxy_http_version 1.1;
        proxy_set_header Upgrade \$http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host \$host;
        proxy_read_timeout 86400;
    }

    # Health check endpoint
    location /health {
        access_log off;
        return 200 "healthy\n";
        add_header Content-Type text/plain;
    }
}
EOF

    # Enable the site
    ln -sf /etc/nginx/sites-available/streamlit /etc/nginx/sites-enabled/
    
    # Remove default site
    rm -f /etc/nginx/sites-enabled/default
    
    # Test Nginx configuration
    print_info "Testing Nginx configuration..."
    nginx -t
    
    # Restart Nginx
    print_info "Restarting Nginx..."
    systemctl restart nginx
    systemctl enable nginx
    
    print_info "Nginx configured successfully!"
    print_info "Application will be accessible at: http://$SERVER_NAME"
}

################################################################################
# Step 4: Set Up SSL Certificate (Optional)
################################################################################

setup_ssl() {
    print_header "Step 4: Setting Up SSL Certificate"
    
    if [ -z "$DOMAIN" ]; then
        print_warning "No domain specified. Skipping SSL setup."
        print_info "To set up SSL later with a domain, run:"
        print_info "  sudo certbot --nginx -d your-domain.com"
        return
    fi
    
    print_info "Setting up SSL certificate for $DOMAIN..."
    print_warning "Make sure your domain DNS is pointing to this server!"
    
    read -p "Continue with SSL setup? (y/n): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Skipping SSL setup."
        return
    fi
    
    # Get SSL certificate
    certbot --nginx -d "$DOMAIN" --non-interactive --agree-tos --register-unsafely-without-email
    
    print_info "SSL certificate installed successfully!"
    print_info "Application will be accessible at: https://$DOMAIN"
}

################################################################################
# Step 5: Create Password Hashing Utility
################################################################################

create_password_utility() {
    print_header "Step 5: Creating Password Hashing Utility"
    
    HASH_SCRIPT="$APP_DIR/hash_password.py"
    
    cat > "$HASH_SCRIPT" << 'EOF'
#!/usr/bin/env python3
"""
Password Hashing Utility for OAM IITK Application

Usage:
    python3 hash_password.py your_password
    
This will output the SHA-256 hash of your password to use in secrets.toml
"""

import sys
import hashlib

def hash_password(password):
    """Hash a password using SHA-256"""
    return hashlib.sha256(password.encode()).hexdigest()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 hash_password.py <password>")
        sys.exit(1)
    
    password = sys.argv[1]
    hashed = hash_password(password)
    
    print(f"\nOriginal password: {password}")
    print(f"Hashed password:   {hashed}")
    print(f"\nAdd this to your ~/.streamlit/secrets.toml:")
    print(f'username = "{hashed}"')
    print()
EOF

    chmod +x "$HASH_SCRIPT"
    
    print_info "Password hashing utility created at: $HASH_SCRIPT"
    print_info "Usage: python3 $HASH_SCRIPT your_password"
}

################################################################################
# Step 6: Create Enhanced Authentication Module
################################################################################

create_auth_module() {
    print_header "Step 6: Creating Enhanced Authentication Module"
    
    AUTH_FILE="$APP_DIR/funcs/stauth_secure.py"
    
    cat > "$AUTH_FILE" << 'EOF'
"""
Enhanced Streamlit Authentication with Security Features
- Password hashing (SHA-256)
- Session timeout
- Login attempt limiting
"""

import streamlit as st
import hashlib
import time

# Configuration
SESSION_TIMEOUT = 3600  # 1 hour in seconds
MAX_LOGIN_ATTEMPTS = 5
LOCKOUT_TIME = 300  # 5 minutes in seconds

def hash_password(password):
    """Hash a password using SHA-256"""
    return hashlib.sha256(password.encode()).hexdigest()

def check_session_timeout():
    """Check if session has timed out"""
    if "last_activity" not in st.session_state:
        st.session_state.last_activity = time.time()
        return False
    
    current_time = time.time()
    if current_time - st.session_state.last_activity > SESSION_TIMEOUT:
        # Session expired
        st.session_state.clear()
        return True
    
    # Update last activity time
    st.session_state.last_activity = current_time
    return False

def check_lockout():
    """Check if user is locked out due to too many failed attempts"""
    if "login_attempts" not in st.session_state:
        st.session_state.login_attempts = 0
        st.session_state.lockout_time = 0
    
    current_time = time.time()
    
    # Check if still in lockout period
    if st.session_state.login_attempts >= MAX_LOGIN_ATTEMPTS:
        if current_time < st.session_state.lockout_time:
            remaining = int(st.session_state.lockout_time - current_time)
            return True, remaining
        else:
            # Lockout period expired, reset attempts
            st.session_state.login_attempts = 0
            st.session_state.lockout_time = 0
    
    return False, 0

def check_password():
    """Returns `True` if the user had the correct password."""

    def password_entered():
        """Checks whether a password entered by the user is correct."""
        # Check lockout status
        is_locked, remaining = check_lockout()
        if is_locked:
            st.error(f"🔒 Too many failed attempts. Try again in {remaining} seconds.")
            return
        
        username = st.session_state["username"]
        password = st.session_state["password"]
        
        # Hash the entered password
        hashed_password = hash_password(password)
        
        # Check against hashed passwords in secrets
        if username in st.secrets["passwords"] and hashed_password == st.secrets["passwords"][username]:
            st.session_state["password_correct"] = True
            st.session_state["authenticated_user"] = username
            st.session_state.login_attempts = 0  # Reset attempts on success
            st.session_state.last_activity = time.time()
            del st.session_state["password"]  # Don't store password
            del st.session_state["username"]  # Don't store username
        else:
            st.session_state["password_correct"] = False
            st.session_state.login_attempts += 1
            
            # Set lockout time if max attempts reached
            if st.session_state.login_attempts >= MAX_LOGIN_ATTEMPTS:
                st.session_state.lockout_time = time.time() + LOCKOUT_TIME

    # Check for session timeout
    if "password_correct" in st.session_state and st.session_state["password_correct"]:
        if check_session_timeout():
            st.warning("⏱️ Session expired. Please log in again.")
            return False

    # Check lockout status
    is_locked, remaining = check_lockout()
    if is_locked:
        st.error(f"🔒 Too many failed login attempts. Please try again in {remaining} seconds.")
        time.sleep(1)
        st.rerun()
        return False

    if "password_correct" not in st.session_state:
        # First run, show inputs for username + password
        st.text_input("Username", on_change=password_entered, key="username")
        st.text_input("Password", type="password", on_change=password_entered, key="password")
        
        # Show remaining attempts
        if "login_attempts" in st.session_state and st.session_state.login_attempts > 0:
            remaining_attempts = MAX_LOGIN_ATTEMPTS - st.session_state.login_attempts
            st.warning(f"⚠️ {remaining_attempts} attempts remaining")
        
        return False
    elif not st.session_state["password_correct"]:
        # Password incorrect, show input + error
        st.text_input("Username", on_change=password_entered, key="username")
        st.text_input("Password", type="password", on_change=password_entered, key="password")
        st.error("😕 User not known or password incorrect")
        
        # Show remaining attempts
        remaining_attempts = MAX_LOGIN_ATTEMPTS - st.session_state.login_attempts
        if remaining_attempts > 0:
            st.warning(f"⚠️ {remaining_attempts} attempts remaining")
        
        return False
    else:
        # Password correct
        return True
EOF

    print_info "Enhanced authentication module created at: $AUTH_FILE"
    print_info "To use it, update your app.py import:"
    print_info "  from funcs.stauth_secure import check_password"
}

################################################################################
# Step 7: Set File Permissions
################################################################################

set_permissions() {
    print_header "Step 7: Setting Secure File Permissions"
    
    # Get the actual user (not root)
    ACTUAL_USER=$(logname 2>/dev/null || echo $SUDO_USER)
    
    if [ -z "$ACTUAL_USER" ]; then
        print_warning "Could not determine actual user. Skipping permission setup."
        return
    fi
    
    STREAMLIT_DIR="/home/$ACTUAL_USER/.streamlit"
    
    if [ -d "$STREAMLIT_DIR" ]; then
        print_info "Setting permissions for $STREAMLIT_DIR..."
        chown -R "$ACTUAL_USER:$ACTUAL_USER" "$STREAMLIT_DIR"
        chmod 700 "$STREAMLIT_DIR"
        
        if [ -f "$STREAMLIT_DIR/secrets.toml" ]; then
            chmod 600 "$STREAMLIT_DIR/secrets.toml"
            print_info "Secured secrets.toml (600)"
        fi
    fi
    
    print_info "File permissions set successfully!"
}

################################################################################
# Step 8: Create Monitoring Script
################################################################################

create_monitoring_script() {
    print_header "Step 8: Creating Monitoring Script"
    
    MONITOR_SCRIPT="$APP_DIR/monitor_app.sh"
    
    cat > "$MONITOR_SCRIPT" << 'EOF'
#!/bin/bash

# Monitor Streamlit Application
# Checks if the app is running and restarts if needed

APP_DIR="$HOME/oam_iitk"
LOG_FILE="$APP_DIR/monitor.log"

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >> "$LOG_FILE"
}

# Check if Streamlit is running
if ! pgrep -f "streamlit run app.py" > /dev/null; then
    log_message "Streamlit is not running. Starting..."
    
    cd "$APP_DIR"
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate oam_env
    
    nohup streamlit run app.py > streamlit.log 2>&1 &
    
    log_message "Streamlit started with PID $!"
else
    log_message "Streamlit is running normally"
fi

# Check if Nginx is running
if ! systemctl is-active --quiet nginx; then
    log_message "Nginx is not running. Starting..."
    sudo systemctl start nginx
    log_message "Nginx started"
fi
EOF

    chmod +x "$MONITOR_SCRIPT"
    
    print_info "Monitoring script created at: $MONITOR_SCRIPT"
    print_info "To set up automatic monitoring, add to crontab:"
    print_info "  */5 * * * * $MONITOR_SCRIPT"
}

################################################################################
# Step 9: Display Summary and Next Steps
################################################################################

display_summary() {
    print_header "Security Setup Complete!"
    
    SERVER_IP=$(curl -s ifconfig.me)
    
    echo -e "${GREEN}✓ Nginx reverse proxy configured${NC}"
    echo -e "${GREEN}✓ Firewall rules set up${NC}"
    echo -e "${GREEN}✓ Security headers added${NC}"
    echo -e "${GREEN}✓ Password hashing utility created${NC}"
    echo -e "${GREEN}✓ Enhanced authentication module created${NC}"
    echo -e "${GREEN}✓ Monitoring script created${NC}"
    
    print_header "Next Steps"
    
    echo "1. Hash your passwords:"
    echo "   python3 $APP_DIR/hash_password.py your_password"
    echo ""
    echo "2. Update ~/.streamlit/secrets.toml with hashed passwords:"
    echo "   [passwords]"
    echo "   admin = \"hashed_password_here\""
    echo ""
    echo "3. Update app.py to use enhanced authentication:"
    echo "   from funcs.stauth_secure import check_password"
    echo ""
    echo "4. Start your Streamlit application:"
    echo "   cd $APP_DIR"
    echo "   conda activate oam_env"
    echo "   streamlit run app.py"
    echo ""
    echo "5. Access your application:"
    if [ -z "$DOMAIN" ]; then
        echo "   http://$SERVER_IP"
    else
        echo "   https://$DOMAIN"
    fi
    echo ""
    echo "6. (Optional) Set up monitoring cron job:"
    echo "   crontab -e"
    echo "   Add: */5 * * * * $APP_DIR/monitor_app.sh"
    echo ""
    
    print_warning "Important Security Notes:"
    echo "- Direct access to port $STREAMLIT_PORT is blocked from outside"
    echo "- All traffic goes through Nginx on ports 80/443"
    echo "- Keep your system updated: sudo apt update && sudo apt upgrade"
    echo "- Regularly review access logs: sudo tail -f /var/log/nginx/access.log"
    echo ""
}

################################################################################
# Main Execution
################################################################################

main() {
    print_header "OAM IITK Security Setup Script"
    
    # Check if running as root
    check_root
    
    # Get domain if user wants SSL
    read -p "Do you have a domain name for SSL? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        read -p "Enter your domain name: " DOMAIN
    fi
    
    # Run setup steps
    install_packages
    configure_firewall
    configure_nginx
    setup_ssl
    create_password_utility
    create_auth_module
    set_permissions
    create_monitoring_script
    display_summary
    
    print_info "Security setup completed successfully!"
}

# Run main function
main
