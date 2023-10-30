function email_radar_status_v1(recipients)
% THIS WILL EMAIL  AODN FILE STATUS
% usage:
% email_radar_status_v1(recipients)
% recipients is cell array of email addersses
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','E_mail','yasha.hetzel@gmail.com');
setpref('Internet','SMTP_Username','yasha.hetzel');
setpref('Internet','SMTP_Password','jcfzrdqbqlpghpaq');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');


email_address=recipients;
attachments={['/home2/yasha/radar/data/FV00_AODN_last-VECTOR-files_' datestr(now,'YYYY') '.txt'],['/home2/yasha/radar/data/FV00_AODN_last-RADIAL-files_' datestr(now,'YYYY') '.txt']};
subject=['Latest AODN file status']

disp('Sending email...')
sendmail(email_address,subject,['Latest AODN file status ' datestr(now) ' email yasha.hetzel@uwa.edu to be removed from list' ],attachments)

for i=1:length(recipients)
disp(['Email sent to: ' recipients{i} ])
end

end

