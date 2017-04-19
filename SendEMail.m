function sucessflag = SendEMail
try
mail='BGUMedicalImagingLab@gmail.com';
password='Tammy''sLab';

% since these must be entered LITERALLY, you could create a P-file and load
% the variables by a call to it...

setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props=java.lang.System.getProperties;
pp=props.setProperty('mail.smtp.auth','true'); %#ok
pp=props.setProperty('mail.smtp.socketFactory.class','javax.net.ssl.SSLSocketFactory'); %#ok
pp=props.setProperty('mail.smtp.socketFactory.port','465'); %#ok
sucessflag = true;
catch err
sucessflag = false;
end

% sendmail('zahi_hershko@yahoo.com','Test email', 'Test');